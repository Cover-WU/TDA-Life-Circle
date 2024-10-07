# Package Import ----------------------------------------------------------
# clear the variables in the momory
rm(list = ls())

# import the R packages
library(TDA)
library(MASS)
library(parallel)
library(rgdal) # have resigned since 2023
library(rgeos)
library(raster)
library(Rsolnp) 
library(lubridate)
library(msr)
library(tidyverse)

# below packages are not called any more
# library(HDoutliers)
# library(magick)
# library(data.tree)

# Variables & Constant --------------------------------------------------

# administrative area of Shenzhen
Shenzhen <- readOGR(dsn='./data/shenzhen_polygon', layer = 'Shenzhen_polygon_new', 
                    use_iconv = TRUE, encoding = "UTF-8")

# coordinate system:
# UTM projection with distance unit in meters
# WGS coordniates with longitudes and latitudes
crs.utm <- CRS("+init=epsg:32650")
crs.wgs <- CRS("+proj=longlat +datum=WGS84")

# bootstrapping times and confidence level
B <- 300; alpha=0.05

# set the surplus
SurplusConstant <- pi / 3 # the unit: cubic kilometers.
DistantTreshold <- 1000 # the unit: meters

# proportion of PH
m.o <- 0.05
# batch size of one bootstrap sampling
BatchSize <- 3000
# spatial resolution at two stage.
ResoUrban <- 500
ResoNeigh <- 50

# whether the environment is linux OS
is.Linux <- TRUE
# only extract the component surrounding the residence
is.Base <- FALSE
# number of parallel processes
nCores <- detectCores()

# origin of research period (for community-level analysis)
timeline.origin <- ymd_hms("2019-10-31 23:00:00")


# Data Loading Area -------------------------------------------------------

# reading data for community life circles
file.path <- './data/quinze_juin_dataset/'
full.list <- list.files(path = file.path, pattern = "*.csv") %>% as.list

read.func <- function(filename){
  full.path <- paste(file.path, filename, sep = '/')
  data <- read_csv(full.path, col_types = cols(t_start=col_datetime(), t_end=col_datetime())) %>% 
    select(c(commune, X, Y, t_start, t_end))
  return(data)
}

location <- map(full.list, read.func)
location.names <- full.list %>% sapply(function(x) gsub('.csv', '', x))

# convert data frame to spatial data frame
convert2SpDf <- function(data){
  coord <- data %>% select(c(X, Y)) %>% as.matrix
  data <- data %>% select(-c(X, Y))
  SpDF <- SpatialPointsDataFrame(coord, data, proj4string = crs.utm)
  return(SpDF)}

names(location) <- location.names
location.Points <- lapply(location, convert2SpDf)
names(location.Points) <- location.names

# convert the coordination system from UTF to WGS
location.Points.Vis <- location.Points %>% lapply(spTransform, crs.wgs)

# Functional Area ---------------------------------------------------------

# Given a data matrix, return the bound box. Padding the width when necessary.
boundBox <- function(datamat, pad = 1000){
  ld.corner <- datamat %>% apply(2, min) %>% {. - pad}
  ur.corner <- datamat %>% apply(2, max) %>% {. + pad}
  bound.box <- rbind(ld.corner, ur.corner)
  row.names(bound.box) <- c("min", "max")
  return(bound.box)
}

# Slicing the dataset into batches with given size. Smaller batch weights less.
batchPartition <- function(datamat, batch.size = BatchSize, random.seed = 5){
  # for repetition
  set.seed(random.seed)
  N <- nrow(datamat)
  # shuffle the dataset
  datamat <- datamat[sample(N),]
  # label the batch sequential number
  batch.code <- ceiling(seq(N) / batch.size)
  # split the dataset into batches
  batches <- split(datamat, batch.code) %>% map(~matrix(., byrow=FALSE, ncol=ncol(datamat)))
  # calculate the smallest weight
  last <- rev(batches)[[1]]
  weight <- nrow(last) / batch.size
  # weight is a number less than one.
  return(list(batches=batches, weight=weight))
}


# find one topological component at one valley in DTM topography.
fluidTopoCup <- function(rasteri, thres, minCellID){
  # subset the area above a distance threshold.
  raster.cut <- calc(rasteri, function(x) ifelse(x > thres, 0, thres - x))
  # count and label the component of the raster.
  raster.component <- clump(raster.cut, directions=4)
  # get the label of target minimum.
  component.label <- raster.component[minCellID]
  # select the component containing the target label.
  raster.mask <- match(raster.component, component.label)
  # mask spatial analysis
  raster.cut.component <- raster::mask(raster.cut, raster.mask)
  # return the raster layer
  return(raster.cut.component)
}

# Whether surplus of a component is smaller than given one.
# unit is cubic kilometer
atVolume <- function(raster.cut.component, Target){
  # compute the approximate integration value.
  integ <- cellStats(raster.cut.component, stat = 'sum', na.rm=TRUE)
  # correct the unit.
  integ.std <- integ / 1000 * prod(res(raster.cut.component)/1000)
  # judge: should the current threshold be smaller, or greater?
  judge <- integ.std <= Target
  return(judge)
}


# create Morse-Smale complex to partition the ascending manifolds
MorseSmaleComplexSeg <- function(rasteri, siglevel){
  # calculate the value
  y <- values(rasteri)
  xcoords <- coordinates(rasteri)
  # compute the Morse-Smale complex
  ms.complex <- msc.nn(-y, xcoords, knn = 8, pLevel = siglevel)
  # todo: add a warning reminder when number of asending manifold is not equal to significant topo signals.
  
  # extract the result and generate segmentation raster layer
  ascending <- ms.complex$level[[1]]$ascending
  raster.seg <- init(rasteri, ascending)
  # return the result
  return(raster.seg)
}


# extract the topological component at given cell ID with threshold designated.
searchCupCityLevel <- function(rasteri, Thres, cellID, absolute=T){
  # find the pixel of the closest value
  minima <- ifelse(absolute, 0, rasteri[cellID])
  # subset the area above a distance threshold.
  thres.cut <- minima + Thres
  # must be relative threashold above local minima
  if(!absolute || rasteri[cellID] < thres.cut){
    res <- fluidTopoCup(rasteri, thres.cut, cellID)
  }
  # else: absolution height
  else{
    res <- raster(rasteri); res[] <- NA
  }
  return(res)
}


# cross sectional area at city level.
# this function differs to the former without any extraction
sectCupCityLevel <- function(rasteri, Thres){
  # subset the area above a distance threshold.
  res <- clamp(rasteri, upper=Thres, useValue=FALSE)
  return(res)
}

# extract the valley cup indicated by Morse-Smale complex.
# rasteri is the original field value raster.
plexCupCityLevel <- function(rasteri, raster.seg, Thres, cellID){
  # extract segmentation label
  seg.tag <- raster.seg[cellID]
  # match the segmentation area
  raster.match <- match(raster.seg, seg.tag)
  # mask the original distance raster
  raster.clip <- mask(rasteri, raster.match)
  # return the search result
  raster.cup <- searchCupCityLevel(raster.clip, Thres, cellID, absolute=F)
  return(raster.cup)
}


# search small cup through interval binary search with certain surplus.
searchOneSmallCup <- function(rasteri, Target, cellID){
  # find the pixel of the closest value
  left.close <- rasteri[cellID]
  right.close <- maxValue(rasteri)
  interval <- c(left.close, right.close)
  # 1k iterations at most
  for(iter in 1:1000){
    if(diff(interval) < 0.1 | iter == 1000){
      result <- fluidTopoCup(rasteri, interval[1], cellID)
      return(result)
    }
    midpoint <- mean(interval)
    if(atVolume(fluidTopoCup(rasteri, midpoint, cellID), Target)) interval[1] <- midpoint else interval[2] <- midpoint
    msgline <- sprintf("Iteration %d, Closed Interval [%.1f, %.1f]...\n", iter, interval[1], interval[2])
    cat(msgline)
  }
}

# search many cups with Lagrange optimization.
# `cellID` is a vector
searchManySmallCups <- function(rasteri, Target, cellID){
  while (TRUE) {
    # initialized theta parameters
    theta <- rasteri[cellID] + 0.01
    # optimized parameters by optimizer
    optimParams <- thresOptimLagrange(rasteri, theta, cellID, Target)
    # return a combined raster with multiple parts.
    compolist <- map2(optimParams, cellID, fluidTopoCup, rasteri=rasteri)
    
    popidx <- vector()
    cat('----------------------------------------------------------\n')
    cat('checking the connectedness...\n')
    for(idx in seq(length(cellID))){
      # traverse the component
      raster.current <- compolist[[idx]]
      # label the spatial range of that component
      myfunc <- function(x){ifelse(is.na(x), 0, 1)}
      r.overlay <- calc(raster.current, fun=myfunc)
      # for the first component: keep it as the template
      if(isTRUE(all.equal(idx, 1))) {template <- r.overlay}
      # for later components: if not overlapped, union it with the current template
      else if(maxValue(r.overlay + template) < 2) {template <- template + r.overlay}
      # if overlapped, record and merge it later.
      else {popidx <- append(popidx, idx); cat(sprintf("Compo at cell %d merged.\n", cellID[idx]))}
    }
    if(!length(popidx)){
      # if no overlapping, overlay all components and return
      circ.life <- compolist %>% reduce(cover)
      return(circ.life)
    }
    else {
      # else, pop the overlapping components and restart the optimization...
      cellID <- cellID[-popidx]
      # except for only one component left.
      if(length(cellID) == 1) return(searchOneSmallCup(rasteri, Target, cellID))
      cat('Cells above are merged, we solve the equation again.\n\n')
    }
  }
}

# Lagrangian optimizer
thresOptimLagrange <- function(rasteri, theta, minCellID, Target){
  # equality constraint: the total amount of liquid is constant.
  equalitySurplusConstraint <- function(theta, rasteri, minCellID){
    # number of components.
    betti <- length(theta)
    # prepare the vector.
    V <- vector('numeric', betti)
    # calculate the surplus of each one.
    for(i in 1:betti){
      theta.i <- theta[i]
      compo <- fluidTopoCup(rasteri, theta.i, minCellID[i])
      V[i] <- cellStats(compo, 'sum', na.rm=TRUE) / 1000 * prod(res(compo)/1000)
    }
    return(sum(V))
  }
  # optimization target: the variation between the cups should be minimized.
  optimSurplusTarget <- function(theta, rasteri, minCellID){
    # number of components.
    betti <- length(theta)
    # prepare the vector.
    V <- vector('numeric', betti)
    # calculate the surplus of each one.
    for(i in 1:betti){
      theta.i <- theta[i]
      compo <- fluidTopoCup(rasteri, theta.i, minCellID[i])
      V[i] <- cellStats(compo, 'sum', na.rm=TRUE) / 1000 * prod(res(compo)/1000)
    }
    return(var(V))
  }
  
  # solve the Lagrange multiplier problem
  solution <- solnp(theta, fun=optimSurplusTarget, eqfun=equalitySurplusConstraint, eqB=Target, rasteri=rasteri, minCellID=minCellID)
  return(solution$pars)
}

# Accelerate the DTM computation for big dataset.
# The strategy is batch tricks.
dtmBig <- function(X, grid, m.o, weight){
  # for smaller dataset, just run original function.
  if(nrow(X) <= BatchSize){
    res <- dtm(X, grid, m.o, weight = weight)
  }
  # or partition random batches.
  else{
    has.weight <- length(weight) > 1
    # if has weight, append the weight vector behind the dataset
    # and partition them together.
    if(has.weight){
      X <- cbind(X, weight)
    }
    batch.partition <- batchPartition(X)
    # weights or no weights, the X have n or n+1 columns.
    X.batches <- batch.partition$batches
    # small batch size weight of the remainder batch, it has no association with point weight at all.
    rem.weight <- batch.partition$weight
    weight.vec <- rep_len(1, length(X.batches))
    weight.vec[length(weight.vec)] <- rem.weight
    if(has.weight){
      # split the X and w to recover.
      w.batches <- lapply(X.batches, function(x) x[, ncol(x)])
      X.batches <- lapply(X.batches, function(x) x[,-ncol(x)])
      # mapply the DTM function
      dtm.result <- mapply(dtm, X = X.batches, weight=w.batches, 
                           MoreArgs = list(Grid=grid, m0=m.o), SIMPLIFY = FALSE)
    }
    else{
      # directly lapply the DTM
      dtm.result <- lapply(X.batches, dtm, Grid=grid, m0=m.o)
    }
    # average by batch size weight.
    dtm.agg <- Reduce(rbind, dtm.result) * weight.vec
    res <- dtm.agg %>% apply(2, sum) %>% {. / sum(weight.vec)}
  }
  return(res)
}

# draw confidence band of the persistent diagram.
# the same logic of batch partition as before.
bootstrapDiagramBig <- function(X, FUN, lim, by, maxdimension = length(lim) / 2 - 1,
                                sublevel = TRUE, library = "GUDHI", B = 30, alpha = 0.05,
                                distance = "bottleneck", dimension = min(1, maxdimension),
                                p = 1, parallel = FALSE, printProgress = FALSE, weight = NULL,
                                isLinux = TRUE, ...){
  # Here we find a bug: when the weight argument is supplemented, 
  # the bootstrapping cannot feedback any valid band (always zero).
  if(nrow(X) <= BatchSize){
    set.seed(5)
    res <- bootstrapDiagram(X, FUN, lim, by, maxdimension, sublevel, 
                            library, B, alpha, distance, dimension,
                            p, parallel, printProgress, weight,
                            ...)
  }
  else{
    has.weight <- length(weight) > 1
    if(has.weight){
      X <- cbind(X, weight)
    }
    batch.X <- batchPartition(X)
    rem.weight <- batch.X$weight
    batch.X <- batch.X$batches
    weight.vec <- rep_len(1, length(batch.X))
    weight.vec[length(weight.vec)] <- rem.weight
    if(has.weight){
      batch.w <- lapply(batch.X, function(x) x[, ncol(x)])
      batch.X <- lapply(batch.X, function(x) x[,-ncol(x)])
      bands <- mapply(bootstrapDiagram, X=batch.X, weight=batch.w, 
                      MoreArgs = list(FUN=FUN, lim=lim, by=by, maxdimension=maxdimension,
                                      sublevel=sublevel, library=library, B=B, alpha=alpha,
                                      distance=distance, dimension=dimension, p=p,
                                      parallel=parallel, printProgress=printProgress, ...=...),
                      SIMPLIFY = TRUE)
    }
    else{
      bands <- sapply(batch.X, bootstrapDiagram, FUN, lim, by, maxdimension, sublevel, 
                      library, B, alpha, distance, dimension, p, parallel, printProgress, 
                      weight, ...)
    }
    res <- weighted.mean(bands, weight.vec)
  }
  return(res)
}

# whether the two cell is in the same connected component.
isCellConnected <- function(r, c0, c1){
  r.clump <- clump(r, directions=4)
  return(r.clump[c0] == r.clump[c1])
}


# search the connecting point of two cells.
searchConnectPoint <- function(r, c0, c1, left.close, right.close){
  interval <- c(left.close, right.close)
  for(iter in 1:1000){
    if(diff(interval) < 1e-3 | iter == 1000){
      return(interval[1])
    }
    midpoint <- mean(interval)
    r.cut <- clamp(r, upper=midpoint, useValue=F)
    if(isCellConnected(r.cut, c0, c1)) interval[2] <- midpoint else interval[1] <- midpoint
    msgline <- sprintf("Iteration %d, Closed Interval [%.3f, %.3f]...\n", iter, interval[1], interval[2])
    cat(msgline)
  }
}

# parallel optimization: zig-zag strategy, 12344321
zigzagSort <- function(x, n){
  x <- x[order(x, decreasing = TRUE)]
  sortvec <- rep(c(seq(1,n),seq(n, 1)), length = length(x))
  sortvec <- order(sortvec)
  x <- x[sortvec]
  return(x)
}

# parallel optimization: stepwise strategy, 
stepSort <- function(x, n){
  x <- sort(x, decreasing = T)
  t <- floor(length(x) / n) 
  sortvec <- rep(c(seq(1,t)), length = length(x))
  sortvec <- order(sortvec)
  x <- x[sortvec]
  return(x)
}

# Community: Detection ----------------------------------------------------


CityLevelLC <- function(location.cloud, m.o, weight=FALSE){
  # delineate the spatial extent of city level life circle.
  # whether time duration should be considered as weight when DTM calculating.
  if(weight){
    timep <- location.cloud %>% transmute(timep = interval(t_start, t_end)) %>% .$timep
    weight.vec <- timep / hours(1)
    weight.boot <- weight.vec
  }
  else{
    weight.vec <- 1
    weight.boot <- NULL
  }
  
  # turn data frame to matrix
  location.cloud <- location.cloud %>% select(c(X, Y)) %>% as.matrix()
  
  # # deprecated:
  # if(noise.screen){
  #   outlier.index <- HDoutliers(location.cloud, alpha = .05)
  #   location.cloud <- location.cloud[outlier.index, ]
  #   if(weight) weight.vec <- weight.vec[outlier.index, ]
  # }
  
  # set up the DTM configuration.
  grid.box <- boundBox(location.cloud)
  lim <- grid.box
  resolution <- ResoUrban
  x.seq <- seq(grid.box["min", "X"], grid.box["max", "X"], by=resolution)
  y.seq <- seq(grid.box["min", "Y"], grid.box["max", "Y"], by=resolution)
  grid <- expand.grid(x.seq, y.seq)
  
  gc()
  # compute DTM
  dtm.grid <- dtmBig(location.cloud, grid, m.o, weight=weight.vec)
  # compute persistent diagram.
  grid.diag <- gridDiag(location.cloud, dtm, lim=lim, by=resolution, library = "Dionysus", location = TRUE, 
                        printProgress=TRUE, m0=m.o, weight=weight.vec)
  # band computed from PD.
  band0 <- bootstrapDiagramBig(location.cloud, dtm, lim, by=resolution, maxdimension = 1, B=B, alpha=alpha, 
                               dimension=0, parallel=T, printProgress = TRUE, m0=m.o, weight=weight.boot, 
                               isLinux = is.Linux)
  gc()
  # extract the stat 
  df.PD <- grid.diag$diagram %>% .[,,drop=F] %>% as_tibble() %>% filter(dimension==0)
  df.PDim <- df.PD %>% mutate(Lifetime = Death-Birth, rn=row_number())
  # filter out the not significant signals.
  df.PD.reduce <- df.PDim %>% filter(Lifetime > 2*band0)
  if(nrow(df.PD.reduce) == 0) df.PD.reduce <- df.PDim %>% slice_max(Lifetime, n=1)
  
  # representative component (topological clusters)
  beta.0.idx <- df.PD.reduce  %>% .$rn
  if(!is.vector(grid.diag$birthLocation)) 
    beta.0.loc <- grid.diag$birthLocation[beta.0.idx, , drop=FALSE]
  else beta.0.loc <- matrix(grid.diag$birthLocation, nrow = 1)
  
  # Convert the matrix to a georeferenced raster.
  z.dtm <- matrix(dtm.grid, ncol=length(y.seq), nrow=length(x.seq))
  raster.dtm <- raster(z.dtm) %>% t() %>% flip()
  extent(raster.dtm) <- c(min(x.seq), max(x.seq), min(y.seq), max(y.seq)) + resolution * c(-1,1,-1,1) / 2
  crs(raster.dtm) <- crs.utm
  
  # location of cluster cores.
  local.coord <- beta.0.loc
  local.cell <- cellFromXY(raster.dtm, local.coord)
  
  # distance thresholding for all cores to delineate the cluster extent.
  raster.life.circle.global <- lapply(local.cell, searchCupCityLevel, Thres=DistantTreshold, rasteri=raster.dtm) 
  raster.life.circle.combined <- raster.life.circle.global %>% reduce(cover)
  return(list(m.o=m.o, PD=grid.diag, band=band0, resRstr=raster.life.circle.combined, oriRstr=raster.dtm, 
              repCompo=local.cell, CompoList=raster.life.circle.global))
}

# higher level wrapper of the fine raster function.
City2Commune <- function(city.res, act.cloud, m.o, base=FALSE, weight=FALSE){
  # extract the resultant raster
  raster.life.circle.combined <- city.res$resRstr
  # find the components across the raster.
  raster.life.compo <- clump(raster.life.circle.combined, direction=4)
  # extract only home component:
  if(base){
    # by default, the deepest component is life circle surrounding the residence and located at the first.
    tick <- raster.life.compo[city.res$repCompo[1]]
    # apply the `Coarse2Fine` to life circle to get finer shape.
    res <- Coarse2Fine(raster.city = raster.life.compo, i=tick, act.cloud=act.cloud, m.o = m.o, weight=weight)
  }
  # or extract all components:
  else{
    # count the components
    number.compo <- cellStats(raster.life.compo, 'max', na.rm=TRUE)
    # apply the `Coarse2Fine` to each compo
    res <- lapply(seq(number.compo), Coarse2Fine, raster.city=raster.life.compo, 
                  act.cloud=act.cloud, m.o=m.o, weight=weight)
  }
  return(res)
}

City2CommuneOnWrite <- function(city.res, act.cloud, m.o, base=FALSE, weight=FALSE, path='.'){
  # writing version
  res <- City2Commune(city.res, act.cloud, m.o, base, weight)
  save(res, file = path)
}

# this function return a fine grained raster of each component domain.
Coarse2Fine <- function(raster.city, i, act.cloud, m.o, weight){
  # i is the serial number of one life component.
  
  # select one component and get its spatial extent
  raster.life.subcompo <- match(raster.city, i)
  raster.life.subcompo <- trim(raster.life.subcompo)
  bound.box <- extent(raster.life.subcompo)
  # convert the coarse extent into a polygon for later subsetting.
  poly.life.subcompo <- rasterToPolygons(raster.life.subcompo) %>% gUnaryUnion()
  
  # find the act points fallen in (subsetting)
  point.location.subcloud <- act.cloud[poly.life.subcompo, ]
  location.subcloud <- point.location.subcloud@coords
  
  if(weight){
    timep <- point.location.subcloud@data %>% transmute(timep = interval(t_start, t_end)) %>% .$timep
    weight.vec <- timep / hours(1)
    weight.boot <- weight.vec
  }
  else{
    weight.vec <- 1
    weight.boot <- NULL
  }
  
  # set fine-grained resolution.
  resolution <- ResoNeigh
  
  # generate the grid: grid point is the center of the pixels.
  # bound box padding
  xpadding <- 2 * xres(raster.city); ypadding <- 2 * yres(raster.city)
  # center offsetting.
  x.seq <- seq(attr(bound.box,"xmin") - xpadding + resolution/2, attr(bound.box,"xmax") + xpadding - resolution/2, by=resolution)
  y.seq <- seq(attr(bound.box,"ymin") - ypadding + resolution/2, attr(bound.box,"ymax") + ypadding - resolution/2, by=resolution)
  grid <- expand.grid(x.seq, y.seq)
  
  # TDA process
  # limit the area.
  bound.box <- bound.box %>% as.matrix(nrow=2) %>% t()
  padding <- c(-1,1) %*% t(c(xpadding, ypadding))
  lim <- bound.box + padding
  # tuning the best m0 ratio param.
  gc()  # recycle
  
  # calculate the distance surface with the best m0.
  dtm.grid <- dtmBig(location.subcloud, grid, m.o, weight=weight.vec)
  # surface to matrix
  z.dtm <- matrix(dtm.grid, ncol=length(y.seq), nrow=length(x.seq))
  # matrix to raster
  raster.commune <- flip(raster(t(z.dtm)))
  # geocode the raster
  extent(raster.commune) <- c(min(x.seq), max(x.seq), min(y.seq), max(y.seq)) + resolution/2 * c(-1, 1,-1, 1)
  crs(raster.commune) <- crs.utm
  
  # generate the persistent diagram
  grid.diag <- gridDiag(location.subcloud, dtm, lim=lim, by=resolution, library = "Dionysus", 
                        location = TRUE, printProgress=TRUE, m0=m.o, weight=weight.vec)
  
  # bands on the zero and first dimension
  band0 <- bootstrapDiagramBig(location.subcloud, dtm, lim, by=resolution, maxdimension = 1, B=B, 
                               alpha=alpha, dimension=0, parallel=T, weight=weight.boot,
                               printProgress = T, m0=m.o, isLinux = is.Linux)
  gc()
  band1 <- bootstrapDiagramBig(location.subcloud, dtm, lim, by=resolution, maxdimension = 1, B=B, 
                               alpha=alpha, dimension=1, parallel=T, weight=weight.boot,
                               printProgress = T, m0=m.o, isLinux = is.Linux)
  gc()
  
  df.PD <- grid.diag$diagram %>% .[,,drop=F] %>% as_tibble()
  if(all(df.PD$dimension == 0)) band1prime <- NULL else band1prime <- band1
  df.PDim <- df.PD %>% mutate(Lifetime = Death-Birth, rn=row_number()) %>% group_by(dimension) %>% nest()
  # filter out the not significant signals.
  df.PD.reduce <- df.PDim %>% ungroup() %>%  mutate(data = map2(data, 2*c(band0, band1prime), ~filter(.x, Lifetime > .y))) %>% unnest(cols=c(data))
  
  # representative component (topological clusters)
  beta.0.idx <- df.PD.reduce %>% filter(dimension==0) %>% .$rn
  if(!is.vector(grid.diag$birthLocation)) 
    beta.0.loc <- grid.diag$birthLocation[beta.0.idx, , drop=FALSE]
  else beta.0.loc <- matrix(grid.diag$birthLocation, nrow = 1)
  cellID <- cellFromXY(raster.commune, beta.0.loc)
  
  # representative loop (cavities)
  beta.1.idx <- df.PD.reduce %>% filter(dimension==1) %>% .$rn
  beta.1.loop <- grid.diag$cycleLocation[beta.1.idx]
  sigLoops <- beta.1.loop
  
  return(list(topoSignal=df.PD.reduce, rasterFine=raster.commune, m.o=m.o, repreCompo=list(cell=cellID, xy=beta.0.loc), 
              repreLoop=sigLoops, band=list(zero=band0, one=band1), gridDiag=grid.diag, 
              coarseShape=poly.life.subcompo, subsetActs=point.location.subcloud))
}

# compute circles in just one area
Commune2Circ <- function(fine.res, Target){
  
  # how many separated components in the area?
  betti.o <- length(fine.res$repreCompo$cell)
  # if there's only one
  if(betti.o == 1){
    circ.search <- searchOneSmallCup(fine.res$rasterFine, Target=Target, cellID=fine.res$repreCompo$cell)
    return(circ.search)
  }
  # if there's more than two
  else{
    circ.search <- searchManySmallCups(fine.res$rasterFine, Target=Target, cellID=fine.res$repreCompo$cell)
    return(circ.search)
  }
}

# Stage 1: city level.
# input the variables into the processes.
cl.export <- c("dtmBig", "bootstrapDiagramBig", "boundBox", "batchPartition", "searchCupCityLevel",
               "fluidTopoCup", "B", "alpha", "ResoUrban", "crs.utm", "BatchSize", "is.Linux", "DistantTreshold")
# parse the parallel processes.
parCluster <- makeCluster(nCore)
# set up packages in each process.
clusterEvalQ(parCluster, {library(TDA); library(raster); library(tidyverse)})
clusterExport(parCluster, cl.export)
# parallel cluster running.
tryCatch(
  {cityLevelLCRes <- parLapplyLB(parCluster, location, CityLevelLC, m.o=m.o)},
  finally = {stopCluster(parCluster)})


# Stage 2: neighborhood level DTM raster computation.
cl.export <- c("dtmBig", "bootstrapDiagramBig", "boundBox", "batchPartition", "Coarse2Fine", "B", "alpha", 
               "ResoNeigh", "crs.utm", "BatchSize", "is.Linux")
parCluster <- makeCluster(nCores)
clusterEvalQ(parCluster, {library(TDA); library(raster); library(tidyverse); 
  library(rgeos)})
clusterExport(parCluster, cl.export)
# 
tryCatch(
  {fineCommuneRes <- clusterMap(parCluster, City2Commune, cityLevelLCRes, location.Points, 
                                MoreArgs = list(m.o=m.o, base=is.Base), .scheduling = "dynamic")},
  finally = {stopCluster(parCluster)})

# if not simply the home circles, flat the hierarchy.
if(!is.Base){
  fineCommuneRes <- unlist(fineCommuneRes, recursive = FALSE)
}

# Stage 3:neighborhood level delineation.
cl.export <- c("searchOneSmallCup", "searchManySmallCups", "fluidTopoCup", "atVolume", "thresOptimLagrange")
parCluster <- makeCluster(nCores)
clusterEvalQ(parCluster, {library(raster); library(tidyverse); 
  library(Rsolnp)})
clusterExport(parCluster, cl.export)
tryCatch(
  {CommuneCirc <- parLapplyLB(parCluster, fineCommuneRes, Commune2Circ, Target=SurplusConstant)},
  finally = {stopCluster(parCluster)})



# Split the Activities by Time --------------------------------------------

# set the duration of the time windows.
time.band.width.hours = 2
group.ticks <- seq(1, 23, 2)
# vector: 23, 1, 3, ..., 21
group.ticks <- c(group.ticks[length(group.ticks)], group.ticks[-length(group.ticks)])
# Chinese Shichen name.
names(group.ticks) <- c("Zi", "Chou", "Yin", "Mao", "Chen", "Si",
                        "Wu", "Wei", "Shen", "You", "Xu", "Hai")

# Extract the index of corresponding tibble
time.ticks <- timeline.origin + hours(seq(0, 24 * 30 + 2, by=time.band.width.hours))
time.bands <- interval(time.ticks[-length(time.ticks)], time.ticks[-1])
time.segs <- tibble(band.interv = time.bands) %>% 
  mutate(timehash = row_number() - 1, 
         grp_level = hour(int_start(band.interv)), 
         grp_label = names(group.ticks)[match(grp_level, group.ticks)])
rm(time.ticks, time.bands)

# split a point dataset into list by time window.
extractIntersect <- function(location.df){
  res <- location.df %>% mutate(act.interv = t_start %--% t_end,
                                key.start = floor((timeline.origin %--% t_start) %/% hours(time.band.width.hours)), 
                                key.end = floor((timeline.origin %--% t_end) %/% hours(time.band.width.hours))) %>% 
    select(-c(t_start, t_end)) %>% 
    mutate(timehash = map2(key.start, key.end, seq)) %>% 
    unnest(timehash) %>% 
    select(-c(key.start, key.end)) %>% 
    inner_join(time.segs, by='timehash') %>% 
    mutate(scatter = intersect(act.interv, band.interv),
           t_start = int_start(scatter), 
           t_end = int_end(scatter)) %>% 
    select(-c(band.interv, act.interv, timehash, scatter)) %>% 
    mutate(grp_label = factor(grp_label, levels = names(group.ticks))) %>% 
    split(.$grp_label) %>% 
    map(~ select(., -grp_label))
  return(res)
}

# apply the split
location.timesplit <- lapply(location, extractIntersect) 
# flat the list hierarchy.
location.timesplit.parallel <- location.timesplit %>% unlist(recursive = FALSE)
# coordnation system conversion.
location.timesplit.parallel.spDf <- lapply(location.timesplit.parallel, convert2SpDf)
# get the problem scale.
subproblem.scale <- sapply(location.timesplit.parallel, nrow)

# Community Temporal Evolution ------------------------------------------------------

# adjust the bootstrap iteration number to fit the data scale.
B <- 300; is.Base <- TRUE
BatchSize <- 3000

# Running the delineation in parallel.
cl.export <- c("dtmBig", "bootstrapDiagramBig", "boundBox", "batchPartition", "searchCupCityLevel",
               "fluidTopoCup", "B", "alpha", "ResoUrban", "crs.utm", "BatchSize", "is.Linux", "DistantTreshold")
parCluster <- makeCluster(nCores)
clusterEvalQ(parCluster, {library(TDA); library(raster); library(tidyverse); library(lubridate)})
clusterExport(parCluster, cl.export)
tryCatch(
  {cityLevelLCResTimeSplit <- parLapplyLB(parCluster, location.timesplit.parallel, CityLevelLC, m.o=m.o, weight=TRUE)},
  finally = {stopCluster(parCluster)})


# change the order to optimize running speed (stepwise strategies)
scale.order <- stepSort(x = subproblem.scale, n = nCores)
cityLevelLCResTimeSplit <- cityLevelLCResTimeSplit[names(scale.order)]
location.timesplit.parallel.spDf <- location.timesplit.parallel.spDf[names(scale.order)]

cl.export <- c("dtmBig", "bootstrapDiagramBig", "boundBox", "batchPartition", "Coarse2Fine", "B", "alpha", 
               "ResoNeigh", "crs.utm", "BatchSize", "is.Linux")
parCluster <- makeCluster(nCores)
clusterEvalQ(parCluster, {library(TDA); library(raster); library(tidyverse); library(lubridate);
  library(rgeos)})
clusterEvalQ(parCluster, sink(paste0("./parallel_log/logging_", Sys.getpid(), ".txt")))
clusterExport(parCluster, cl.export)
tryCatch(
  {fineCommuneResTimeSplit <- clusterMap(parCluster, City2Commune, cityLevelLCResTimeSplit, location.timesplit.parallel.spDf, 
                                         MoreArgs = list(m.o=m.o, base=is.Base, weight=TRUE), .scheduling = "dynamic")},
  finally = {stopCluster(parCluster); save(fineCommuneResTimeSplit, file='./fineCommuneRes_05_TimeSplit_b3k.RData')})


cl.export <- c("searchOneSmallCup", "searchManySmallCups", "fluidTopoCup", "atVolume", "thresOptimLagrange")
parCluster <- makeCluster(nCores)
clusterEvalQ(parCluster, {library(raster); library(tidyverse); library(Rsolnp)})
clusterEvalQ(parCluster, sink(paste0("./optimization_log/logging_", Sys.getpid(), ".txt")))
clusterExport(parCluster, cl.export)
tryCatch(
  {CommuneCircTimeSplit <- parLapplyLB(parCluster, fineCommuneResTimeSplit, Commune2Circ, Target=SurplusConstant)},
  finally = {stopCluster(parCluster); save(CommuneCircTimeSplit, file='./CommuneCirc_05_TimeSplit_b3k.RData')})


# Individual: Inspection --------------------------------------------------

# read individual-level dataset.
file.path <- './data/espace_juin_dataset/'
full.list <- list.files(path = file.path, pattern = "*.csv") %>% as.list

read.func <- function(filename){
  full.path <- paste(file.path, filename, sep = '/')
  data <- read_csv(full.path, col_types = cols(t_start=col_datetime(), t_end=col_datetime())) %>% 
    select(c(who, X, Y, t_start, t_end))
  return(data)
}

# variable names are repeated as before, so that the content will be covered.
location <- map(full.list, read.func)
location.names <- full.list %>% sapply(function(x) {x <- rev(unlist(strsplit(x, '_', fixed=T)))[1]; gsub('.csv', '', x)})

names(location) <- location.names
location.Points <- lapply(location, convert2SpDf)
names(location.Points) <- location.names

# note: the difference between `Vis` and `vIs` for community-level and individual-level respectively.
location.Points.vIs <- location.Points %>% lapply(spTransform, crs.wgs)

# not simply the home circle but all anchor points across the city.
is.Base <- FALSE
Shenzhen.utm <- spTransform(Shenzhen, CRSobj = crs.utm)

# get the problem scale.
subproblem.scale <- sapply(location, nrow)
subproblem.scale


# rewrite the city level life circle `CityLevelLC` code, but modify small points:
# - change the grid box spanning
# - introduce the Morse-Smale complex and ascending manifolds to divide the range of the anchor 
CityLevelAnchor <- function(location.cloud, m.o, weight=FALSE){
  if(weight){
    timep <- location.cloud %>% transmute(timep = interval(t_start, t_end)) %>% .$timep
    weight.vec <- timep / hours(1)
    weight.boot <- weight.vec
  }
  else{
    weight.vec <- 1
    weight.boot <- NULL
  }
  
  location.cloud <- location.cloud %>% select(c(X, Y)) %>% as.matrix()
  
  # use the bound box of the whole Shenzhen area.
  grid.box <- bbox(Shenzhen.utm) %>% t()
  colnames(grid.box) <- colnames(grid.box) %>% toupper()
  lim <- grid.box
  resolution <- ResoUrban
  x.seq <- seq(grid.box["min", "X"], grid.box["max", "X"], by=resolution)
  y.seq <- seq(grid.box["min", "Y"], grid.box["max", "Y"], by=resolution)
  grid <- expand.grid(x.seq, y.seq)
  
  gc()
  dtm.grid <- dtmBig(location.cloud, grid, m.o, weight=weight.vec)
  grid.diag <- gridDiag(location.cloud, dtm, lim=lim, by=resolution, library = "Dionysus", location = TRUE, 
                        printProgress=TRUE, m0=m.o, weight=weight.vec)

  band0 <- bootstrapDiagramBig(location.cloud, dtm, lim, by=resolution, maxdimension = 1, B=B, alpha=alpha, 
                               dimension=0, parallel=F, printProgress = T, m0=m.o, weight=weight.boot, 
                               isLinux = is.Linux)
  gc()
  df.PD <- grid.diag$diagram %>% .[,,drop=F] %>% as_tibble() %>% filter(dimension==0)
  df.PDim <- df.PD %>% mutate(Lifetime = Death-Birth, rn=row_number())
  df.PD.reduce <- df.PDim %>% filter(Lifetime > 2*band0)
  if(nrow(df.PD.reduce) == 0) df.PD.reduce <- df.PDim %>% slice_max(Lifetime, n=1)
  
  beta.0.idx <- df.PD.reduce  %>% .$rn
  if(!is.vector(grid.diag$birthLocation)) 
    beta.0.loc <- grid.diag$birthLocation[beta.0.idx, , drop=FALSE]
  else beta.0.loc <- matrix(grid.diag$birthLocation, nrow = 1)
  
  z.dtm <- matrix(dtm.grid, ncol=length(y.seq), nrow=length(x.seq))
  raster.dtm <- raster(z.dtm) %>% t() %>% flip()
  extent(raster.dtm) <- c(min(x.seq), max(x.seq), min(y.seq), max(y.seq)) + resolution * c(-1,1,-1,1) / 2
  crs(raster.dtm) <- crs.utm
  
  local.coord <- beta.0.loc
  local.cell <- cellFromXY(raster.dtm, local.coord)
  
  # Morse-Smale complex computation.
  raster.seg <- MorseSmaleComplexSeg(raster.dtm, siglevel=2*band0)
  raster.life.circle.global <- lapply(local.cell, plexCupCityLevel, Thres=DistantTreshold, rasteri=raster.dtm, raster.seg=raster.seg) 
  raster.life.circle.combined <- raster.life.circle.global %>% reduce(cover)
  return(list(m.o=m.o, PD=grid.diag, band=band0, resRstr=raster.life.circle.combined, oriRstr=raster.dtm, 
              segRstr=raster.seg, repCompo=local.cell, CompoList=raster.life.circle.global))
}


# wrapper of the Anchor function: the number of the nearest neighbors rather than the ratio is fixed.
CityLevelAS <- function(location.cloud, k, weight=FALSE){
  # use eps to prevent in wrong integer rounding.
  eps <- 1e-7
  m0 <- k / nrow(location.cloud) - eps
  # m0 should be unity at most.
  m0 <- ifelse(m0 > 1, 1, m0)
  res <-  CityLevelAnchor(location.cloud, m.o = m0, weight = weight)
  return(res)
}

# modified from `Coarse2Fine`, with fixed number of nearest neighbors,
# besides, the params `raster.city` becomes a list of raster rather than the
# one clump raster
CoarseAnc2FineStruct <- function(raster.city, i, act.cloud, k, weight){
  # i is the serial number of one anchor.
  
  # select one component and get its spatial extent
  raster.life.subcompo <- raster.city[[i]]
  raster.life.subcompo <- trim(raster.life.subcompo)
  bound.box <- extent(raster.life.subcompo)
  poly.life.subcompo <- rasterToPolygons(raster.life.subcompo) %>% gUnaryUnion()
  
  # do not clip the data, here we use all activity records
  point.location.subcloud <- act.cloud
  location.subcloud <- point.location.subcloud@coords

  if(weight){
    timep <- point.location.subcloud@data %>% transmute(timep = interval(t_start, t_end)) %>% .$timep
    weight.vec <- timep / hours(1)
    weight.boot <- weight.vec
  }
  else{
    weight.vec <- 1  
    weight.boot <- NULL
  }
  
  xrange <- attr(bound.box, 'xmax') - attr(bound.box, 'xmin')
  yrange <- attr(bound.box, 'ymax') - attr(bound.box, 'ymin')
  xdim <- xrange / ResoNeigh
  ydim <- yrange / ResoNeigh
  
  # if the component is too large, we externally set the upper limit of the raster length
  # in order to accelarate the computation.
  max.dim <- 200
  if(max(xdim, ydim) > max.dim){
    resolution <- ceiling(max(xrange, yrange) / max.dim)
  }
  else{
    resolution <- ResoNeigh
  }

  # given k asceratain the m0
  eps <- 1e-7
  m0 <- k / nrow(location.subcloud) - eps
  m0 <- ifelse(m0 > 1, 1, m0)
  m.o <- m0
  
  # grid cetting
  xpadding <- 2 * xres(raster.life.subcompo); ypadding <- 2 * yres(raster.life.subcompo)
  x.seq <- seq(attr(bound.box,"xmin") - xpadding + resolution/2, attr(bound.box,"xmax") + xpadding - resolution/2, by=resolution)
  y.seq <- seq(attr(bound.box,"ymin") - ypadding + resolution/2, attr(bound.box,"ymax") + ypadding - resolution/2, by=resolution)
  grid <- expand.grid(x.seq, y.seq)
  
  # TDA process
  bound.box <- bound.box %>% as.matrix(nrow=2) %>% t()
  padding <- c(-1,1) %*% t(c(xpadding, ypadding))
  lim <- bound.box + padding
  gc()  # recycle
  
  dtm.grid <- dtmBig(location.subcloud, grid, m.o, weight=weight.vec)
  z.dtm <- matrix(dtm.grid, ncol=length(y.seq), nrow=length(x.seq))
  raster.commune <- flip(raster(t(z.dtm)))
  extent(raster.commune) <- c(min(x.seq), max(x.seq), min(y.seq), max(y.seq)) + resolution/2 * c(-1, 1,-1, 1)
  crs(raster.commune) <- crs.utm
  
  # generate the PD and its bands.
  grid.diag <- gridDiag(location.subcloud, dtm, lim=lim, by=resolution, library = "Dionysus", 
                        location = TRUE, printProgress=TRUE, m0=m.o, weight=weight.vec)
  
  band0 <- bootstrapDiagramBig(location.subcloud, dtm, lim, by=resolution, maxdimension = 1, B=B, 
                               alpha=alpha, dimension=0, parallel=F, weight=weight.boot,
                               printProgress = T, m0=m.o, isLinux = is.Linux)
  gc()
  band1 <- bootstrapDiagramBig(location.subcloud, dtm, lim, by=resolution, maxdimension = 1, B=B, 
                               alpha=alpha, dimension=1, parallel=F, weight=weight.boot,
                               printProgress = T, m0=m.o, isLinux = is.Linux)
  gc()
  
  df.PD <- grid.diag$diagram %>% .[,,drop=F] %>% as_tibble()
  
  if(all(df.PD$dimension == 0)) band1prime <- NULL else band1prime <- band1
  df.PDim <- df.PD %>% mutate(Lifetime = Death-Birth, rn=row_number()) %>% group_by(dimension) %>% nest()
  df.PD.reduce <- df.PDim %>% ungroup() %>%  mutate(data = map2(data, 2*c(band0, band1prime), ~filter(.x, Lifetime > .y))) %>% unnest(cols=c(data))
  if(nrow(df.PD.reduce) == 0) df.PD.reduce <- df.PDim %>% filter(dimension==0) %>% slice_max(Lifetime, n=1)
  
  # representative component
  beta.0.idx <- df.PD.reduce %>% filter(dimension==0) %>% .$rn
  if(!is.vector(grid.diag$birthLocation)) 
    beta.0.loc <- grid.diag$birthLocation[beta.0.idx, , drop=FALSE]
  else beta.0.loc <- matrix(grid.diag$birthLocation, nrow = 1)
  cellID <- cellFromXY(raster.commune, beta.0.loc)
  
  # representative loop
  beta.1.idx <- df.PD.reduce %>% filter(dimension==1) %>% .$rn
  beta.1.loop <- grid.diag$cycleLocation[beta.1.idx]
  sigLoops <- beta.1.loop
  
  # clip the raster by the polygon.
  raster.commune <- mask(raster.commune, poly.life.subcompo)
  valid.topo.idx <- which(!is.na(raster.commune[cellID]))
  cellID <- cellID[valid.topo.idx]
  beta.0.loc <- beta.0.loc[valid.topo.idx, ]
  
  return(list(topoSignal=df.PD.reduce, rasterFine=raster.commune, m.o=m.o, repreCompo=list(cell=cellID, xy=beta.0.loc), 
              repreLoop=sigLoops, band=list(zero=band0, one=band1), gridDiag=grid.diag, 
              coarseShape=poly.life.subcompo, subsetActs=point.location.subcloud))
}

# high functional wrapper of the anchor finer delineation.
Close2Anchor <- function(city.res, act.cloud, k, base=FALSE, weight=FALSE){
  # high resolution DTM raster within coarse extent
  raster.clipper <- city.res$resRstr
  # morse-smale complex segmentation raster
  raster.segmented <- city.res$segRstr
  # high resolution original DTM raster
  raster.source <- city.res$oriRstr
  # clip the segmentation by extent
  raster.clip <- mask(raster.segmented, raster.clipper)
  # unique value of segmentation (number of valleys)
  compo.unique <- unique(raster.clip)
  raster.compo.list <- lapply(compo.unique, function(x){
    # match the anchor with the number
    anc.match <- match(raster.clip, x)
    # return to extract the raster
    res.sub <- mask(raster.source, anc.match)
    return(res.sub)
  })
  number.compo <- length(raster.compo.list)
  res <- lapply(seq(number.compo), CoarseAnc2FineStruct, raster.city=raster.compo.list,
                act.cloud=act.cloud, k=k, weight=weight)
  
  return(res)
}


# compute circles in just one area
Commune2Circ <- function(fine.res, Target){
  # how many separated components in the area?
  betti.o <- length(fine.res$repreCompo$cell)
  # then we multiply the number to render each component with nearly the same target value.
  Target <- Target * betti.o
  # if there's only one
  if(betti.o == 1){
    circ.search <- searchOneSmallCup(fine.res$rasterFine, Target=Target, cellID=fine.res$repreCompo$cell)
    return(circ.search)
  }
  # if there's more than two
  else{
    circ.search <- searchManySmallCups(fine.res$rasterFine, Target=Target, cellID=fine.res$repreCompo$cell)
    return(circ.search)
  }
}

# wrap the function for parallel error
safeWrapper <- function(func, errorPrefix = "Error:") {
  safeFunc <- function(...) {
    tryCatch({
      result <- do.call(func, args = list(...))
      return(result)
    }, error = function(e) {
      errorMsg <- paste(errorPrefix, conditionMessage(e))
      return(errorMsg)
    })
  }
  return(safeFunc)
}

# set the confidence level to unity (high sensitivity)
alpha <- 1
# Running the delineation in parallel.
cl.export <- c("dtmBig", "bootstrapDiagramBig", "boundBox", "batchPartition", "MorseSmaleComplexSeg", 
               "plexCupCityLevel", "safeWrapper", 
               "searchCupCityLevel", "CityLevelAnchor","CityLevelAS", "fluidTopoCup", "B", "alpha", 
               "ResoUrban", "crs.utm", "BatchSize", "is.Linux", "DistantTreshold", "Shenzhen.utm")
parCluster <- makeCluster(nCores)
clusterEvalQ(parCluster, sink(paste0("./parallel_log/logging_", Sys.getpid(), ".txt")))
clusterEvalQ(parCluster, {library(TDA); library(raster); library(msr); library(tidyverse)})
clusterExport(parCluster, cl.export)
tryCatch(
  {cityLevelAncRes <- parLapplyLB(parCluster, location, safeWrapper(CityLevelAS), k=10)},
  finally = {stopCluster(parCluster); save(cityLevelAncRes, file='./individualCityAnchor_k10.RData')})


cl.export <- c("dtmBig", "bootstrapDiagramBig", "boundBox", "batchPartition", "CoarseAnc2FineStruct", "B", "alpha",
               "ResoNeigh", "crs.utm", "BatchSize", "is.Linux", "Close2Anchor")
parCluster <- makeCluster(nCores)
clusterEvalQ(parCluster, sink(paste0("./parallel_log/logging_", Sys.getpid(), ".txt")))
clusterEvalQ(parCluster, {library(TDA); library(raster); library(tidyverse);
  library(rgeos)})
clusterExport(parCluster, cl.export)
tryCatch(
  {fineAnchorRes <- clusterMap(parCluster, safeWrapper(Close2Anchor), cityLevelAncRes, location.Points,
                               MoreArgs = list(k=10, base=is.Base), .scheduling = 'dynamic')},
  finally = {stopCluster(parCluster); save(fineAnchorRes, file='./individualCommunityAnchor_k10.RData')})


if(!is.Base){
  fineCommuneRes <- unlist(fineAnchorRes, recursive = FALSE)
}

# rename the results
left <- names(fineCommuneRes) %>% substr(1, 9)
right <- names(fineCommuneRes) %>% substring(10)
names(fineCommuneRes) <- paste(left, right, sep = '.')


subproblem.scale <- sapply(fineCommuneRes, function(x) length(x$repreCompo$cell))
scale.order <- stepSort(x = subproblem.scale, n = nCores)
fineCommuneRes <- fineCommuneRes[names(scale.order)]


cl.export <- c("searchOneSmallCup", "searchManySmallCups", "fluidTopoCup", "atVolume", "thresOptimLagrange",
               "safeWrapper", "Commune2Circ")
parCluster <- makeCluster(nCores)
clusterEvalQ(parCluster, sink(paste0("./parallel_log/logging_", Sys.getpid(), ".txt")))
clusterEvalQ(parCluster, {library(raster); library(tidyverse); 
  library(Rsolnp)})
clusterExport(parCluster, cl.export)
tryCatch(
  {CommuneCirc <- parLapplyLB(parCluster, fineCommuneRes, safeWrapper(Commune2Circ), Target=SurplusConstant)},
  finally = {save(CommuneCirc, file='./individualStructAnchor.RData'); stopCluster(parCluster)})



# Visualization -----------------------------------------------------------

# Points of activities descriptive visualization
library(ggspatial)
library(ggmap)
library(ggrgl)
library(devoutrgl)
library(scico)
library(png)
library(ggnewscale)
library(RColorBrewer)
library(animation)
library(aspace)
library(rayshader)
library(sf)

# canvas bound box
shenzhen.bbox <- c(113.7, 22.4, 114.7, 22.9)

# get the basemap
# set the stadia map api key
stadia_apikey <- "b4972941-a168-4aa0-b341-7b34c776c3e7"
register_stadiamaps(stadia_apikey) 
sz.basemap <- get_stadiamap(bbox = shenzhen.bbox, zoom = 11, maptype = 'stamen_terrain')

# draw basemap in ggmap
base.image <- ggmap(sz.basemap, darken = 0.6) + 
  # add north arrow
  annotation_north_arrow(location='tr', style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl") +
  layer_spatial(Shenzhen, fill='white', alpha=0.3) +
  theme_void()
# output the map into png image.
base.path <- './product/basemap.png'
ggsave(base.path, width = 10, height = 5.43, dpi=300)
base.raster <- readPNG(base.path)

# the polygon of communties.
Communes35 <- readOGR(dsn='./data/communities', layer = 'communes35', 
                      use_iconv = TRUE, encoding = "UTF-8")

# function to convert UTM coordnates to WGS.
point2wgs <- function(x) {xy_transform(x[,1], x[,2], from = 32650)}

# colormap of the two analysis
color.vec.commune <- rainbow(35)
color.vec.individ <- rainbow(20, s=0.5, v=0.5)

# create long data frame for community 3D visualization, without internal geocoordinates.
sorge.sub <- function(df.sub) {
  attribute <- df.sub@data
  loccoords <- df.sub@coords
  return(cbind(attribute, loccoords))
}

# Community activity dataset.
act.segments <- location.Points.Vis %>% lapply(sorge.sub) %>% reduce(rbind)
act.segments <- act.segments %>% mutate(z_start = 60 * hour(t_start) + minute(t_start) + round(second(t_start) / 60),
              z_end = 60 * hour(t_end) + minute(t_start) + round(second(t_start) / 60),
              commune = if_else(str_starts(commune, "桃源村"), "桃源村", commune))


# set the 3D box axis (or box edge)
tridbox <- shenzhen.bbox + 0 * c(1,1,-1,-1)
box3d <- tibble(
  x = rep(tridbox[c(1, 3)], 4),
  xend = c(rep(tridbox[c(1, 3)], 2), rep(tridbox[c(3, 1)],each=2) ),
  y = rep(tridbox[c(2, 4, 2, 4)], each=2),
  yend =c(rep(tridbox[c(2, 4)],each=2), rep(tridbox[c(2, 4)], 2)),
  z = rep(c(0, 1440), each=4),
  zend = rep(1440, 8),
)

# We also tried 3D rectangle but failed:
# geom_rect_z(aes(xmin = shenzhen.bbox[1], ymin = shenzhen.bbox[2], xmax=shenzhen.bbox[3], ymax = shenzhen.bbox[4], z=1440), 
#             extrude = T, extrude_face_fill = NA, extrude_edge_colour='black', color='black', fill=NA,
#             material = list(alpha=0), 
#             show.legend = F) 


# collective activity visualization in 3D space
# the 3D engine is ggrgl
g <- ggplot() + 
  annotation_raster(base.raster, xmin = shenzhen.bbox[1],
                                  ymin = shenzhen.bbox[2], xmax=shenzhen.bbox[3],
                                  ymax = shenzhen.bbox[4]) +
  coord_fixed() +
  geom_segment_3d(data=act.segments, aes(x=X, y=Y, z=z_start, xend=X, yend=Y, zend=z_end, color=commune), alpha = 0.2) +
  scale_color_manual(guide = NULL, values = color.vec.commune) +
  geom_segment_3d(data = box3d, aes(x=x, y=y, z=z, xend=xend, yend=yend, zend=zend),  color='black', show.legend = F) +
  scale_x_continuous(name = "Longitude", breaks = seq(shenzhen.bbox[1], shenzhen.bbox[3], 0.2), 
                     limits = shenzhen.bbox[c(1,3)], expand = c(0,0)) +
  scale_y_continuous(name = "Latitude", breaks = seq(shenzhen.bbox[2], shenzhen.bbox[4], 0.1),
                     limits = shenzhen.bbox[c(2,4)], expand = c(0,0)) +
  theme_ggrgl() 

# set the rgl apparatus.
rgldev(filename = './product/describe_acts_activity.png', dpi=300,
       fov = 30, view='flat', view_angle = 30, zscale = 3, zoom=0.8, close_window = T)
g
invisible(dev.off())


# individual level activity visualization in 3D space
act.segments <- location.Points.vIs %>% lapply(sorge.sub) %>% reduce(rbind)
act.segments <- act.segments %>% mutate(z_start = 60 * hour(t_start) + minute(t_start) + round(second(t_start) / 60),
                                        z_end = 60 * hour(t_end) + minute(t_start) + round(second(t_start) / 60),
                                        who = as_factor(who))

g <- ggplot() + 
  annotation_raster(base.raster, xmin = shenzhen.bbox[1],
                    ymin = shenzhen.bbox[2], xmax=shenzhen.bbox[3],
                    ymax = shenzhen.bbox[4]) +
  coord_fixed() +
  geom_segment_3d(data=act.segments, aes(x=X, y=Y, z=z_start, xend=X, yend=Y, zend=z_end, color=who), alpha = 0.2) +
  scale_color_manual(guide = NULL, values = color.vec.individ) +
  geom_segment_3d(data = box3d, aes(x=x, y=y, z=z, xend=xend, yend=yend, zend=zend),  color='black', show.legend = F) +
  scale_x_continuous(name = "Longitude", breaks = seq(shenzhen.bbox[1], shenzhen.bbox[3], 0.2), 
                     limits = shenzhen.bbox[c(1,3)], expand = c(0,0)) +
  scale_y_continuous(name = "Latitude", breaks = seq(shenzhen.bbox[2], shenzhen.bbox[4], 0.1),
                     limits = shenzhen.bbox[c(2,4)], expand = c(0,0)) +
  theme_ggrgl() 

rgldev(filename = './product/describe_acts_anchor.png', dpi=300,
       fov = 30, view='flat', view_angle = 30, zscale = 3, zoom=0.8, close_window = T)
g

invisible(dev.off())

# Activity scatter plot: throughout the November.
# draw activity event in cyan scatter points.
act.scatter.plotter <- function(act.df){
  g <- ggmap(sz.basemap, alpha=0.9, darken = 0.6) +
    layer_spatial(Shenzhen, alpha=0.3) +
    layer_spatial(act.df, color='#8DEEEE', size=0.5, alpha=0.05) +
    annotation_north_arrow(location='tr', style = north_arrow_fancy_orienteering) +
    annotation_scale(location = "bl") +
    xlab(NULL) + ylab(NULL)
  return(g)
}
location.draw <- lapply(location.Points.Vis, act.scatter.plotter)

# Next part we draw the activity space

# At first is the community level life circle
load('./image_result/new_result/CommuneCirc_05.RData')
load('./image_result/new_result/fineCommuneRes_05.RData')

# Draw life circle at fine scale: zoom-in visualization
# extract the community polygon
community.poly <- substr(names(fineCommuneRes), 1, 2) %>% as.integer() %>% 
  lapply(function(d) Communes35[Communes35@data$Encoding == d, ])

# draw community life circle at ith community.
communityTopo <- function(i){
  
  fineRes <- fineCommuneRes[[i]]
  Commune <- CommuneCirc[[i]]
  Reside <- community.poly[[i]]
  
  x <- t(as.matrix(bbox(Commune)))
  # bound box to longtitude and latitude
  bbox.vect <- point2wgs(x) %>% as.matrix() %>% t() %>% c()
  
  compo.loc <- fineRes$repreCompo$xy
  loop <- fineRes$repreLoop
  dfCompo <- point2wgs(compo.loc)
  # draw the base map surrounding the community
  lc.basemap <- get_map(location = bbox.vect, zoom = 13, source = 'stadia', maptype = 'stamen_terrain')
  
  # draw the coarse extent and stay events falling within
  # besides, community life circle and shape
  # cluster cores are pointed out in red triangles
  g <- ggmap(lc.basemap, darken = 0.7) +
    layer_spatial(fineRes$coarseShape %>% spTransform(crs.wgs),
                  color=NA, fill='white', alpha=0.15) +
    layer_spatial(fineRes$subsetActs %>% spTransform(crs.wgs),
                  color='#8DEEEE', size=0.2, alpha=0.1) +
    layer_spatial(Commune, alpha = 0.3) +
    layer_spatial(Reside, color="#444444", fill=NA, linetype='dashed') +
    scale_fill_scico(palette = "lajolla", na.value=NA) +
    geom_point(aes(x=x, y=y), color='red', shape=17, size=2, data=dfCompo)
  
  # if there is a hole within the life circle, highlight it with blue dashed line.
  if(length(loop)) {
    dfLoop <- lapply(loop, function(arr){data.frame(x=arr[,1,1], y=arr[,1,2], xend=arr[,2,1], yend=arr[,2,2])}) %>% reduce(rbind)
    dfLoop <- cbind(point2wgs(dfLoop[,c(1,2)]), point2wgs(dfLoop[,c(3,4)]))
    colnames(dfLoop) <- c('x', 'y', 'xend', 'yend')
    g <- g + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), color="blue", size=1.5, 
                          linejoin='round', linetype='dotted', data = dfLoop)
  }  
  
  # other attachments and decorattions.
  g <- g + annotation_north_arrow(location='tr', style = north_arrow_fancy_orienteering) +
    annotation_scale(location = "bl") +
    xlab(NULL) + ylab(NULL)+
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    guides(fill="none")
  return(g)
}

# lapply the above function and save the figures
lc.draw <- lapply(seq_along(fineCommuneRes), communityTopo)
names(lc.draw) <- names(fineCommuneRes)
save(lc.draw, file='./image_result/topores_community_stat.RData')

# calculate the bounding box in longitudes and latitudes
coordsboxes <- lapply(CommuneCirc, function(x){point2wgs(t(as.matrix(bbox(x)))) %>% as.matrix() %>% t() %>% c()})
# draw life circle figures in Bagualing Dorm and Yuling Garden.
lc.basemap <- get_map(location = coordsboxes$`05_BagualingDorm`, zoom = 13, source = 'stadia', maptype = 'stamen_terrain')
g <- ggmap(lc.basemap, darken = 0.7) +
  xlab(NULL) + ylab(NULL)+
  theme(axis.text.x = element_text(angle=45, hjust=1)) 
ggsave('./product/eg_basemap_bagua.png', width = 7, height = 7, dpi = 300)
plot(g)
lc.basemap <- get_map(location = coordsboxes$`06_YulingGarden`, zoom = 13, source = 'stadia', maptype = 'stamen_terrain')
g <- ggmap(lc.basemap, darken = 0.7) +
  xlab(NULL) + ylab(NULL)+
  theme(axis.text.x = element_text(angle=45, hjust=1))
plot(g)
ggsave('./product/eg_basemap_yuling.png', width = 7, height = 7, dpi = 300)

# draw life circles on all communities.
for(commune.nom in names(lc.draw)){
  plot(lc.draw[[commune.nom]])
  ggsave(paste0('./product/annals_visualization/Communes_LifeCircle_nova/', commune.nom, '.png'), 
         width = 7, height = 7, dpi=300)
  dev.off()
}

# below is drawing all life circles on one single map
# integrate all circles onto the same map.
g <- ggmap(sz.basemap, alpha=0.9, darken = 0.6) +
  layer_spatial(Shenzhen, alpha=0.3) +
  annotation_north_arrow(location='tr', style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl") +
  xlab(NULL) + ylab(NULL)
# create the color map (because one community have two city-level circle)
ans <- c(color.vec.commune[1:3], color.vec.commune[3], color.vec.commune[-(1:3)])
for(i in seq_along(CommuneCirc)){
  # iterate and overlay new circle to basemap with new color.
  g <- g + new_scale_fill() + layer_spatial(CommuneCirc[[i]]) + 
    scale_fill_gradient(low=adjustcolor(ans[i], alpha.f = 0.05), high=ans[i], na.value = NA, guide = 'none')
}
# plot and save
plot(g)
ggsave(paste0('./product/annals_visualization/Communes_LifeCircle_Integ.png'), 
              width = 10, height = 5.45, dpi=600)

# the same as above but with white edge
# yet above circles have transparent edge.
g <- ggmap(sz.basemap, alpha=0.9, darken = 0.6) +
  layer_spatial(Shenzhen, alpha=0.3) +
  annotation_north_arrow(location='tr', style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl") +
  xlab(NULL) + ylab(NULL)
ans <- c(color.vec.commune[1:3], color.vec.commune[3], color.vec.commune[-(1:3)])
for(i in seq_along(CommuneCirc)){
  g <- g + new_scale_fill() + layer_spatial(CommuneCirc[[i]]) + 
    scale_fill_gradient(low='white', high=ans[i], na.value = NA, guide = 'none')
}
plot(g)
ggsave(paste0('./product/annals_visualization/Communes_LifeCircle_Integ_boundary.png'), 
       width = 10, height = 5.45, dpi=600)
dev.off()

# Visualize the communtiy temporal evolution
# load time slices data.
load('./image_result/new_result/fineCommuneRes_05_TimeSplit_b3k.RData')
load('./image_result/new_result/CommuneCirc_05_TimeSplit_b3k.RData')
# Chinese Shichen label
shichen <- c("Zi", "Chou", "Yin", "Mao", "Chen", "Si",
             "Wu", "Wei", "Shen", "You", "Xu", "Hai")
# Decartes product and create names
sorted.name <- expand.grid(shichen, names(location.Points.Vis)) %>% apply(1, function(x) paste(x[2], x[1], sep = '.'))
fineCommuneResTimeSplit <- fineCommuneResTimeSplit[sorted.name]
CommuneCircTimeSplit <- CommuneCircTimeSplit[sorted.name]

# Nested the time slice list.
fineTimeNested <- lapply(names(location.Points.Vis), function(x) fineCommuneResTimeSplit[grepl(x, names(fineCommuneResTimeSplit))])
names(fineTimeNested) <- names(location.Points.Vis)

circTimeNested <- lapply(names(location.Points.Vis), function(x) CommuneCircTimeSplit[grepl(x, names(CommuneCircTimeSplit))])
names(circTimeNested) <- names(location.Points.Vis)

actTimeNested <- lapply(names(location.Points.Vis), function(x) location.timesplit.parallel.spDf[grepl(x, names(location.timesplit.parallel.spDf))])
names(actTimeNested) <- names(location.Points.Vis)

# for each community, query the maximal bound box across the community.
basemapQuery <- function(timeStack){
  boundbox.expand <- rbind(
    map(timeStack, ~ bbox(.$coarseShape) %>% .[,'min']) %>% reduce(rbind) %>% `-`(2 * ResoUrban) %>% apply(2, min),
    map(timeStack, ~ bbox(.$coarseShape) %>% .[,'max']) %>% reduce(rbind) %>% `+`(2 * ResoUrban) %>% apply(2, max)
  )
  basemap.bbox <- point2wgs(boundbox.expand) %>% as.matrix() %>% t() %>% c()
  bm <- get_map(location = basemap.bbox, zoom = 13, source = 'stadia', maptype = 'stamen_terrain')
  return(bm)
}
basemap.bundle <- lapply(fineTimeNested, basemapQuery)
# save the temporal evolution basemap
save(basemap.bundle, file='./image_result/CommunityLifeCircle')
load(file='./image_result/CommunityLifeCircle')

# again, community polygons.
community.poly <- substr(names(location.Points), 1, 2) %>% as.integer() %>% 
  lapply(function(d) Communes35[Communes35@data$Encoding == d, ])
names(community.poly) <- names(location.Points.Vis)

# draw results of a time slice as one frame of the animation
drawOneFrame <- function(fine.slice, circ.slice, act.slice, timeperiod.chr, basemap, reside){
  compo.loc <- fine.slice$repreCompo$xy
  loop <- fine.slice$repreLoop
  dfCompo <- point2wgs(compo.loc)
  
  # set the title
  title <- sprintf("Community Life Circle at %s", timeperiod.chr)
  
  # nearly the same as before
  g <- ggmap(basemap, darken = 0.5) +
    layer_spatial(act.slice %>% spTransform(crs.wgs),
                  color='#8DEEEE', size=0.2, alpha=0.1) +
    layer_spatial(circ.slice, alpha = 0.3) +
    layer_spatial(reside, fill='white', color=NA, alpha=0.3) +
    scale_fill_scico(palette = "lajolla", na.value=NA) +
    geom_point(aes(x=x, y=y), color='red', shape=17, size=2, data=dfCompo)
  
  if(length(loop)) {
    dfLoop <- lapply(loop, function(arr){data.frame(x=arr[,1,1], y=arr[,1,2], xend=arr[,2,1], yend=arr[,2,2])}) %>% reduce(rbind)
    dfLoop <- cbind(point2wgs(dfLoop[,c(1,2)]), point2wgs(dfLoop[,c(3,4)]))
    colnames(dfLoop) <- c('x', 'y', 'xend', 'yend')
    g <- g + geom_segment(aes(x=x, y=y, xend=xend, yend=yend), color="blue", size=1.5, 
                          linejoin='round', linetype='dotted', data = dfLoop)}  
  
  g <- g + annotation_north_arrow(location='tr', style = north_arrow_fancy_orienteering) +
    annotation_scale(location = "bl") +
    xlab(NULL) + ylab(NULL)+
    ggtitle(title) +
    theme(axis.text = element_text(size=rel(1)),
          axis.text.x = element_text(angle=45, hjust=1),
          plot.title = element_text(size=rel(1.75), hjust=0.5)) +
    guides(fill="none")
  return(g)
}


drawAnimation <- function(basemap, reside, fine, circ, act) {
  timeperiods <- c(
    "23:00-1:00", "1:00-3:00", "3:00-5:00", "5:00-7:00",
    "7:00-9:00", "9:00-11:00", "11:00-13:00", "13:00-15:00",
    "15:00-17:00", "17:00-19:00", "19:00-21:00", "21:00-23:00"
  ) 
  # `fine`, `circ`, `act` and `timeperiods` all have length of 12.
  frames <- mapply(drawOneFrame, fine, circ, act, timeperiods, MoreArgs = list(basemap=basemap, reside=reside), SIMPLIFY=FALSE)
  return(frames)
}

# select some exemplary communities when we want to draw the static figures.
figure.communes <- names(basemap.bundle)[c(6, 8, 14)]
basemap.bundle <- basemap.bundle[figure.communes]
community.poly <- community.poly[figure.communes]
fineTimeNested <- fineTimeNested[figure.communes]
circTimeNested <- circTimeNested[figure.communes]
actTimeNested <- actTimeNested[figure.communes]
ans <- mapply(drawAnimation, basemap.bundle, community.poly, fineTimeNested, circTimeNested, actTimeNested, SIMPLIFY=FALSE)

# when drawing pseudo-dynamic maps: export the maps as png images for further visualization.
for(i in seq_along(ans)){
  community.name <- names(ans)[i]
  save.path <- paste0('./product/dyanmic_plots/', community.name, '/')
  dir.create(save.path, recursive = T)
  for(j in seq_along(ans[[community.name]])){
    plot(ans[[community.name]][[j]])
    file.name <- names(ans[[community.name]])[j]
    ggsave(paste0(save.path, sprintf("%02d-", 2*j), file.name, '.png'), height=4, width = 4, dpi=300)
  }
}

# else when exporting the animations.
dev.off()
setwd('./product/annals_visualization/Communes_LifeCircle_Tempo/')
for(i in seq_along(ans)){
  community.name <- names(ans)[i]
  glist <- ans[[community.name]]
  movie.name <- sprintf('%s.gif', community.name)
  saveGIF({for(g in glist) plot(g)}, movie.name=movie.name, interval=1, ani.height=800,
          ani.width= 800)
}
setwd('..'); setwd('..'); setwd('..')

# Visualiza anchors of all individuals
load('./image_result/new_result/individualStructAnchors.RData')
load('./image_result/new_result/individualCommunityAnchor_k10_with_band.RData')

individual <- sort(names(fineAnchorRes))
fineAnchorResFlat <- unlist(fineAnchorRes, recursive = FALSE)
# separate the result name with `.` separator
left <- names(fineAnchorResFlat) %>% substr(1, 9)
right <- names(fineAnchorResFlat) %>% substring(10)
names(fineAnchorResFlat) <- paste(left, right, sep = '.')

# sort the results.
fineAnchorResFlat <- fineAnchorResFlat[sort(names(fineAnchorResFlat))]
IndividAnchor <- CommuneCirc[sort(names(CommuneCirc))]

# draw basemap
gbase <- ggmap(sz.basemap, alpha=0.9, darken = 0.6) +
  layer_spatial(Shenzhen, alpha=0.3) +
  annotation_north_arrow(location='tr', style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl") +
  xlab(NULL) + ylab(NULL)

# draw anchor of individual, whose identifier is named as `ind` arguments.
drawAnchors <- function(ind){
  fine.ind <- fineAnchorResFlat[startsWith(names(fineAnchorResFlat), ind)]
  anchor.ind <- IndividAnchor[startsWith(names(IndividAnchor), ind)]
  acts.ind <- location.Points.vIs[[ind]]
  color.ind <- color.vec.individ[which(individual == ind)]
  # draw the stay events of the `ind` user.
  g <- gbase +  layer_spatial(acts.ind, color='white', size=0.2, alpha=0.1) 
  # iterate and draw all anchors on the map.
  for(i in seq_along(fine.ind)){
    fine.slice <- fine.ind[[i]]
    anchor.pad <- anchor.ind[[i]]
    compo.loc <- fine.slice$repreCompo$xy %>% matrix(ncol=2)
    dfCompo <- point2wgs(compo.loc)
    # fill the anchor with gradual color from white to user colormap
    # label the anchor core with the black triangle.
    g <- g + new_scale_fill() + layer_spatial(anchor.pad) +
      scale_fill_gradient(low = '#ffffff', high = color.ind, 
                          na.value = NA, guide = 'none') + 
      geom_point(aes(x=x, y=y), color='black', shape=17, size=1, data=dfCompo) 
  }
  return(g)
}

# save the visualization
anchor.draw <- lapply(individual, drawAnchors)
names(anchor.draw) <- individual
save(anchor.draw, file='./image_result/topores_anchors.RData')
# print them to png images.
for(nom in names(anchor.draw)){
  plot(anchor.draw[[nom]])
  ggsave(paste0('./product/annals_visualization/Indiv_Anchor_nova/', nom, '.png'), 
         width = 10, height = 5.45, dpi=300)
  dev.off()
}


# Drawing illustration figure for explaining filtration.

# Explaining filtration: generate simulation dataset.
# simulative generating function
gauss <- function(size, mu.x, mu.y, sigma.x, sigma.y=sigma.x, theta=0){
  g <- mvrnorm(n = size, mu=c(0, 0), Sigma = diag(c(sigma.x,sigma.y) ^ 2))
  g <- g %*% matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), ncol = 2) %>% sweep(2, c(mu.x,mu.y), '+')
  return(g)
}

circle <- function(size, x, y, r){
  rho <- sqrt(runif(size))
  theta <- runif(size, max=2*pi)
  jitter <- cbind(r * rho * cos(theta), r * rho * sin(theta)) %>% sweep(2, c(x,y), '+')
  return(jitter)
}

ellipse <- function(size, x, y, a, b, theta=0, jr=0){
  ellip <- circle(size, 0, 0, r=1) %*% diag(c(a, b)) %*% matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), ncol = 2) %>% 
    sweep(2, c(x, y), "+")
  rho <- sqrt(runif(size))
  theta <- runif(size, max=2*pi)
  jitter <- cbind(jr * rho * cos(theta), jr * rho * sin(theta))
  return(ellip + jitter)
}

# bound box area for the university town.
coord.min <- c(113.957252, 22.583541)
coord.max <- c(113.980201, 22.599533)
box.coord <- rbind(coord.min, coord.max)
box.meter <- xy_transform(box.coord[,1], box.coord[,2], to = 32650)

set.seed(5)
simuset <- rbind(
  gauss(300, 113.972971, 22.594547, sigma.x = .0016, sigma.y = .0007, theta = -pi/9), 
  gauss(300, 113.962824, 22.594123, sigma.x = .0018, sigma.y = .0009, theta = 0), 
  gauss(600, 113.966307, 22.589223, sigma.x = .0030, sigma.y = .0008, theta = -pi/9),
  gauss(90, 113.976190, 22.589756, sigma.x = .0008, sigma.y = .0012, theta = 0),
  cbind(runif(60,coord.min[1],coord.max[1]), runif(60, coord.min[2], coord.max[2]))
)
# generation of timestamp: sampling from real dataset.
timedata <- act.segments %>% select(z_start, z_end) %>% sample_n(nrow(simuset), replace = F)
# create spatial dataframe for the simulation dataset.
act.sim <- SpatialPointsDataFrame(coords = simuset, data = timedata, proj4string = crs.wgs)
act.sim.meter <- spTransform(act.sim, CRSobj = crs.utm)
simuset.meter <- act.sim.meter@coords
act.segments <- act.sim.meter@data %>% mutate(X = simuset[,1], Y = simuset[,2])

# basemap for university town.
univ.basemap <- get_stadiamap(bbox = c(coord.min, coord.max), zoom = 15, maptype = 'stamen_terrain')
sim.basemap <- ggmap(univ.basemap) + 
  theme_void()
sim.basemap
# save the basemap as rdata and png.
save(sim.basemap, file='./image_result/simulation_basemap.RData')
save.path <- './product/university_town.png'
ggsave(save.path, width=5, height = 3.78, dpi = 300)
sim.baseimage <- readPNG(save.path)

# similarly draw the 3D edge bound box.
tridbox <- c(coord.min, coord.max) + 0 * c(1,1,-1,-1)
box3d <- tibble(
  x = rep(tridbox[c(1, 3)], 4),
  xend = c(rep(tridbox[c(1, 3)], 2), rep(tridbox[c(3, 1)],each=2) ),
  y = rep(tridbox[c(2, 4, 2, 4)], each=2),
  yend =c(rep(tridbox[c(2, 4)],each=2), rep(tridbox[c(2, 4)], 2)),
  z = rep(c(0, 1440), each=4),
  zend = rep(1440, 8),
)

# draw the 3D segment plot.
g <- ggplot() +
  annotation_raster(sim.baseimage, xmin = coord.min[1], ymin=coord.min[2],
                    xmax = coord.max[1], ymax = coord.max[2]) + 
  geom_segment_3d(data=act.segments, aes(x=X, y=Y, z=z_start, xend=X, yend=Y, zend=z_end), alpha=0.3, color='#ff0000', size=1) +
  geom_segment_3d(data = box3d, aes(x=x, y=y, z=z, xend=xend, yend=yend, zend=zend), color='black', show.legend = F) +
  theme_void()

# adjust the perspective by hand and save the parameters.
# rotation <- rgl::par3d("userMatrix")
# save(rotation, file = './image_result/rotation_3d_Simulation.RData')

# load the perspective settings.
load('./image_result/rotation_3d_Simulation.RData')
# receive viewpoint matrix
rgldev(filename = './product/topoSimulation.png', dpi=300,
       fov = 30, view3d_args = list(userMatrix=rotation), zscale = 3, zoom=1, close_window = T)
g
# close
dev.off()

# draw simulating plots for illustration.
tdaOnSim <- function(simuset.meter, m.o = 0.05){
  alpha <- 0.1
  resolution <- 20
  x.seq <- seq(box.meter[1, 'x'] + resolution/2, box.meter[2, 'x'] - resolution/2, by=resolution)
  y.seq <- seq(box.meter[1, 'y'] + resolution/2, box.meter[2, 'y'] - resolution/2, by=resolution)
  grid <- expand.grid(x.seq, y.seq)
  weight.vec <- 1
  weight.boot <- NULL
  
  lim <- as.matrix(box.meter)
  gc()  # recycle
  
  # calculate the DTM surface
  dtm.grid <- dtmBig(simuset.meter, grid, m.o, weight = weight.vec)
  # surface to matrix
  z.dtm <- matrix(dtm.grid, ncol=length(y.seq), nrow=length(x.seq))
  # matrix to raster
  raster.commune <- flip(raster(t(z.dtm)))
  # georeference the raster
  extent(raster.commune) <- c(min(x.seq), max(x.seq), min(y.seq), max(y.seq)) + resolution/2 * c(-1, 1,-1, 1)
  crs(raster.commune) <- crs.utm

  # generate the persistent diagram with bands
  grid.diag <- gridDiag(simuset.meter, dtm, lim=lim, by=resolution, library = "Dionysus", 
                        location = TRUE, printProgress=TRUE, m0=m.o)
  
  band0 <- bootstrapDiagramBig(simuset.meter, dtm, lim, by=resolution, maxdimension = 1, B=B, 
                               alpha=alpha, dimension=0, parallel=T, weight = weight.boot,
                               printProgress = T, m0=m.o, isLinux = is.Linux)
  gc()
  band1 <- bootstrapDiagramBig(simuset.meter, dtm, lim, by=resolution, maxdimension = 1, B=B, 
                               alpha=alpha, dimension=1, parallel=T, weight = weight.boot,
                               printProgress = T, m0=m.o, isLinux = is.Linux)
  gc()
  
  # filter out the not significant signals.
  df.PD <- grid.diag$diagram %>% .[,,drop=F] %>% as_tibble()
  if(all(df.PD$dimension == 0)) band1prime <- NULL else band1prime <- band1
  df.PDim <- df.PD %>% mutate(Lifetime = Death-Birth, rn=row_number()) %>% group_by(dimension) %>% nest()
  df.PD.reduce <- df.PDim %>% ungroup() %>%  mutate(data = map2(data, 2*c(band0, band1prime), ~filter(.x, Lifetime > .y))) %>% unnest(cols=c(data))
  
  # representative component
  beta.0.idx <- df.PD.reduce %>% filter(dimension==0) %>% .$rn
  if(!is.vector(grid.diag$birthLocation)) 
    beta.0.loc <- grid.diag$birthLocation[beta.0.idx, , drop=FALSE]
  else beta.0.loc <- matrix(grid.diag$birthLocation, nrow = 1)
  cellID <- cellFromXY(raster.commune, beta.0.loc)
  
  # representative loop
  beta.1.idx <- df.PD.reduce %>% filter(dimension==1) %>% .$rn
  beta.1.loop <- grid.diag$cycleLocation[beta.1.idx]
  sigLoops <- beta.1.loop
  
  return(list(topoSignal=df.PD.reduce, rasterFine=raster.commune, m.o=m.o, repreCompo=list(cell=cellID, xy=beta.0.loc), 
              repreLoop=sigLoops, band=list(zero=band0, one=band1), gridDiag=grid.diag))
}

# conduct simulative TDA analysis
tda.sim <- tdaOnSim(simuset.meter)
save(tda.sim, file='./image_result/TDA-simulation.RData')

# draw 3D graph with rayshader package.
elmat <- raster_to_matrix(tda.sim$rasterFine)
filtrate.vert <- c(tda.sim$topoSignal$Birth, tda.sim$topoSignal$Death) %>% sort()
filtrate.vert <- filtrate.vert + 5
shadow <- ray_shade(elmat, zscale = 50, lambert = FALSE)
amb <- ambient_shade(elmat, zscale = 10)

# pseudo-3d visualization of the DTM raster
elmat %>%
  height_shade() %>% 
  add_shadow(shadow, 0.5) %>%
  add_shadow(amb, 0) %>%
  plot_3d(elmat, zscale = 50, fov = 30, theta = -40, phi = 60, 
          windowsize = c(1000, 800), zoom = 0.9)
# save the starting plot
render_snapshot(filename = paste('./product/annals_visualization/Simulate/filtration', 0, '.png'), clear=TRUE)

# visualization of filtration process and save the pictures.
for(i in seq_along(filtrate.vert)){
  elmat %>%
    height_shade() %>% 
    add_shadow(shadow, 0.5) %>%
    add_shadow(amb, 0) %>%
    plot_3d(elmat, zscale = 50, fov = 30, theta = -40, phi = 60, 
            windowsize = c(1000, 800), zoom = 0.9,
            water = TRUE, waterdepth = filtrate.vert[i], wateralpha = 0.5, watercolor = "lightblue",
            waterlinecolor = "white", waterlinealpha = 1)
  render_snapshot(filename = paste('./product/annals_visualization/Simulate/filtration', i, '.png'), clear=TRUE)
}


# Persistent Diagram by ggplot
df.PD <- tda.sim$gridDiag$diagram %>% as.matrix() %>% .[,] %>% as.data.frame()
band0 <- tda.sim$band$zero
band1 <- tda.sim$band$one
gg.PD <- df.PD %>% mutate(Lifetime = Death-Birth, 
                          isSignificant = if_else(dimension==0 & Lifetime > 2*band0 | dimension==1 & Lifetime > 2*band1, T, F),
                          grp = case_when(dimension == 0 & isSignificant ~ 'red',
                                          dimension == 1 & isSignificant ~ 'blue', 
                                          TRUE ~ 'grey'))
stroke.size <- ifelse(gg.PD$grp == 'blue', 2, 0)
limit <- 900
bander <- data.frame(x=c(0,0,limit,limit-band0*2),y=c(band0*2, 0, limit, limit))
ggplot(data = gg.PD) +
  geom_abline(intercept = 0, slope = 1, color='#696969') +
  geom_polygon(aes(x=x, y=y, color=NULL, fill='red', alpha=0.3), data=bander)+
  geom_abline(intercept = 2*band1, slope = 1, linetype='dashed', size=0.8, color='blue')+
  geom_point(aes(x=Birth, y=Death, shape=factor(dimension), size=grp, color=grp, stroke=stroke.size)) +
  scale_color_manual(values = c('#6495ED', '#696969', '#F08080')) + 
  scale_size_manual(values = c(4 ,2, 3))+
  scale_shape_manual(values = c(17, 1, 1))+
  scale_x_continuous(name='Birth', limits = c(0,limit), expand = c(0,0))+
  scale_y_continuous(name='Death', limits = c(0,limit), expand = c(0,0))+
  coord_fixed()+
  theme_bw()+
  theme(axis.title = element_text(face='bold', family = "Myriad"),
        legend.position = 'None')
ggsave('./product/simulation_persistent_diagram.png', width = 6, height = 6, dpi=300)

# Below is simple illustrative plot to generate simple ballistic dataset... 
X <- gauss(1000, 0, 0, 1)
x.seq <- seq(-4, 4, by=0.04)
y.seq <- seq(-4, 4, by=0.04)
grid <- expand.grid(x.seq, y.seq)
dtm.grid <- dtmBig(X, grid, m.o, weight = 1)
z.dtm <- matrix(dtm.grid, ncol=length(y.seq), nrow=length(x.seq))
r.dtm <- flip(raster(t(z.dtm)))
z.dtm <- raster_to_matrix(r.dtm)

shadow <- ray_shade(z.dtm, zscale = 50, lambert = FALSE)
amb <- ambient_shade(z.dtm, zscale = 50)

# explain surplus threshoding...
z.dtm %>% height_shade() %>% 
  add_shadow(shadow, 0.5) %>%
  add_shadow(amb, 0) %>%
  plot_3d(z.dtm, zscale = 0.05, fov = 30, theta = -30, phi = 60, 
          windowsize = c(1000, 800), zoom = 0.9,
          water = TRUE, waterdepth = 0.6, wateralpha = 0.5, watercolor = "lightblue",
          waterlinecolor = "white", waterlinealpha = 1)
render_snapshot(filename = './product/annals_visualization/Simulate/theta_surplus.png', clear=TRUE)

# ... and also distance threshoding
z.dtm %>% height_shade() %>% 
  add_overlay(generate_contour_overlay(z.dtm, levels=0.6, nlevels=1, resolution_multiply = 3, 
                                       color = 'lightblue', linewidth = 10))  %>%
  add_shadow(shadow, 0.5) %>%
  add_shadow(amb, 0) %>%
  plot_3d(z.dtm, zscale = 0.05, fov = 30, theta = -30, phi = 60, 
          windowsize = c(1000, 800), zoom = 0.9,
          water = TRUE, waterdepth = 0.6, wateralpha = 0.1, watercolor = "black")
render_snapshot(filename = './product/annals_visualization/Simulate/theta_dist.png', clear=TRUE)



# The last part: comparison with other measurement.
# Firstly we define a template raster, with fine granularity.
resolution <- 50
box <- round(extent(Shenzhen.utm) / resolution) * resolution
raster.fine.template <- raster(box, crs = crs.utm,
                               nrows = (box@ymax - box@ymin) / resolution,
                               ncols = (box@xmax - box@xmin) / resolution)

# rastarize the point layer from vector data to raster
location.Point <- location.Points.Vis %>% lapply(spTransform, CRSobj=crs.utm)
# choose Bagualing comminty life circle as example case.
act.cloud <- location.Point$`05_BagualingDorm`
act.cloud <- SpatialPoints(act.cloud@coords, proj4string = crs.utm)
raster.fine.count <- rasterize(act.cloud, raster.fine.template, fun='count')

# Kernel density estimation (KDE)
# set the bandwidth of the kernel
kernel.band <- 1000
gauss.kernel <- focalWeight(raster.fine.count, d=kernel.band, type='Gauss')
raster.kde <- focal(raster.fine.count, gauss.kernel, fun=sum, na.rm=T)

# pick out the threshold to filter noise: 0.05
threshold <- 0.05
raster.kde.cut <- calc(raster.kde, function(x) {thres <- threshold; ifelse(x > thres, x - thres, NA)})
# ...and reverse the raster into the polygon spatial layer
kde.delin <- rasterToPolygons(x = raster.kde, n=16, fun = function(x) { x >= threshold })
life.circ.kde <- gBuffer(gUnaryUnion(kde.delin), width=100)

# ggplot figures:
g <- ggmap(sz.basemap, darken = 0.6) +
  layer_spatial(Shenzhen, alpha=0.3) +
  layer_spatial(act.cloud, color='#8DEEEE', size=0.2, alpha=0.01) +
  layer_spatial(raster.kde.cut, alpha=0.5) +
  scale_fill_scico(palette = "lapaz", na.value=NA) +
  layer_spatial(life.circ.kde %>% spTransform(CRS("+proj=longlat +datum=WGS84")), fill=NA, color = 'white') +
  annotation_north_arrow(location='tr', style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl") +
  xlab(NULL) + ylab(NULL) + guides(fill="none")

plot(g)
ggsave('./product/comparison_community/comparison_kde_1km.png', width = 10, height=5.45)

ans <- scico(30, begin = 0, end = 1, palette = "lapaz")
lapaz.blue <- ans[1]

# a finer bandwidth
kernel.band <- 250
gauss.kernel <- focalWeight(raster.fine.count, d=kernel.band, type='Gauss')
raster.kde <- focal(raster.fine.count, gauss.kernel, fun=sum, na.rm=T)

# the threshold holds the same
threshold <- 0.05
raster.kde.cut <- calc(raster.kde, function(x) {thres <- threshold; ifelse(x > thres, x - thres, NA)})
# ...and reverse the raster into the polygon spatial layer
kde.delin <- rasterToPolygons(x = raster.kde, n=16, fun = function(x) { x >= threshold })
# a smaller bandwidth doesn't permit buffering.
life.circ.kde <- gUnaryUnion(kde.delin)

# ggplot figures:
g <- ggmap(sz.basemap, darken = 0.6) +
  layer_spatial(Shenzhen, alpha=0.3) +
  layer_spatial(act.cloud, color='#8DEEEE', size=0.2, alpha=0.01) +
  layer_spatial(raster.kde.cut, alpha=0.5) +
  scale_fill_scico(palette = "lapaz", na.value=NA) +
  layer_spatial(life.circ.kde %>% spTransform(CRS("+proj=longlat +datum=WGS84")), fill=NA, color = 'white') +
  annotation_north_arrow(location='tr', style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl") +
  xlab(NULL) + ylab(NULL) + guides(fill="none")

plot(g)
ggsave('./product/comparison_community/comparison_kde_250m.png', width = 10, height=5.45)

# increase the threshold to 0.5
threshold <- 0.5
raster.kde.cut <- calc(raster.kde, function(x) {thres <- threshold; ifelse(x > thres, x - thres, NA)})
kde.delin <- rasterToPolygons(x = raster.kde, n=16, fun = function(x) { x >= threshold })
life.circ.kde <- gUnaryUnion(kde.delin)

# ggplot figures:
g <- ggmap(sz.basemap, darken = 0.6) +
  layer_spatial(Shenzhen, alpha=0.3) +
  layer_spatial(act.cloud, color='#8DEEEE', size=0.2, alpha=0.01) +
  layer_spatial(raster.kde.cut, alpha=0.5) +
  scale_fill_scico(palette = "lapaz", na.value=NA) +
  layer_spatial(life.circ.kde %>% spTransform(CRS("+proj=longlat +datum=WGS84")), fill=NA, color = 'white') +
  annotation_north_arrow(location='tr', style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl") +
  xlab(NULL) + ylab(NULL) + guides(fill="none")

plot(g)
ggsave('./product/comparison_community/comparison_kde_250m_core.png', width = 10, height=5.45)

# Next appraoch: Standard deviation ellipse (SDE)
ellipse <- calc_sde(points = act.cloud@coords)
ellipse.poly <- st_polygon(list(ellipse$LOCATIONS[, 2:3] %>% as.matrix()))
ellipso.poly.sf <- st_sfc(ellipse.poly, crs = "EPSG:32650")
ellipse.sf <- st_sf(geometry=ellipso.poly.sf)
# ggplot figures:
g <- ggmap(sz.basemap, darken = 0.6) +
  layer_spatial(Shenzhen, alpha=0.3) +
  layer_spatial(act.cloud, color='#8DEEEE', size=0.2, alpha=0.01) +
  layer_spatial(ellipse.sf, color='white', fill=lapaz.blue, alpha=0.5) +
  annotation_north_arrow(location='tr', style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl") +
  xlab(NULL) + ylab(NULL) + guides(fill="none")

plot(g)
ggsave('./product/comparison_community/comparison_ellipse.png', width = 10, height=5.45)

# The next: Minimal convex polygon/hull (MCP)
ch <- chull(act.cloud@coords)
cv <- act.cloud@coords[c(ch, ch[1]), ]
cvxhull <- SpatialPolygons(list(Polygons(list(Polygon(cv)), ID=1)), proj4string = crs.utm)
g <- ggmap(sz.basemap, darken = 0.6) +
  layer_spatial(Shenzhen, alpha=0.3) +
  layer_spatial(act.cloud, color='#8DEEEE', size=0.2, alpha=0.01) +
  layer_spatial(cvxhull, color='white', fill=lapaz.blue, alpha=0.5) +
  annotation_north_arrow(location='tr', style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl") +
  xlab(NULL) + ylab(NULL) + guides(fill="none")

plot(g)
ggsave('./product/comparison_community/comparison_convexhull.png', width = 10, height=5.45)

# Our TDA measurement.
load('./image_result/new_result/CommuneCirc_05.RData')
toporaster <- CommuneCirc$`05`
den.delin <- rasterToPolygons(x = toporaster, n=16)
life.circ.den <- gBuffer(gUnaryUnion(den.delin), width=100)

g <- ggmap(sz.basemap, darken = 0.6) +
  layer_spatial(Shenzhen, alpha=0.3) +
  layer_spatial(act.cloud, color='#8DEEEE', size=0.2, alpha=0.01) +
  layer_spatial(toporaster, alpha=0.8) +
  scale_fill_scico(palette = "lapaz", na.value=NA) +
  layer_spatial(life.circ.den, color='white', fill=NA) +
  annotation_north_arrow(location='tr', style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl") +
  xlab(NULL) + ylab(NULL) + guides(fill="none")

plot(g)
ggsave('./product/comparison_community/comparison_topo.png', width = 10, height=5.45)



# Then is the individual-level comparison
# rastarize the point layer: vector to raster
location.Point <- location.Points.vIs %>% lapply(spTransform, CRSobj=crs.utm)
# We choose one survey user.
act.cloud <- location.Point$`082455786`
act.cloud <- SpatialPoints(act.cloud@coords, proj4string = crs.utm)
raster.fine.count <- rasterize(act.cloud, raster.fine.template, fun='count')

# kernel bandwidth initially set to 1 km
kernel.band <- 1000
gauss.kernel <- focalWeight(raster.fine.count, d=kernel.band, type='Gauss')
raster.kde <- focal(raster.fine.count, gauss.kernel, fun=sum, na.rm=T)

# set the threshold to 0.01
threshold <- 0.01
raster.kde.cut <- calc(raster.kde, function(x) {thres <- threshold; ifelse(x > thres, x - thres, NA)})
# reverse the density raster into the polygon spatial layer
kde.delin <- rasterToPolygons(x = raster.kde, n=16, fun = function(x) { x >= threshold })
life.circ.kde <- gBuffer(gUnaryUnion(kde.delin), width=100)

# ggplot figures:
g <- ggmap(sz.basemap, darken = 0.6) +
  layer_spatial(Shenzhen, alpha=0.3) +
  layer_spatial(act.cloud, color='#8DEEEE', size=0.2, alpha=0.01) +
  layer_spatial(raster.kde.cut, alpha=0.5) +
  scale_fill_scico(palette = "lapaz", na.value=NA) +
  layer_spatial(life.circ.kde %>% spTransform(CRS("+proj=longlat +datum=WGS84")), fill=NA, color = 'white') +
  annotation_north_arrow(location='tr', style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl") +
  xlab(NULL) + ylab(NULL) + guides(fill="none")

plot(g)
ggsave('./product/comparison_individual/comparison_kde_1km.png', width = 10, height=5.45)


# a finer bandwidth
kernel.band <- 250
gauss.kernel <- focalWeight(raster.fine.count, d=kernel.band, type='Gauss')
raster.kde <- focal(raster.fine.count, gauss.kernel, fun=sum, na.rm=T)

# the threshold holds the same
raster.kde.cut <- calc(raster.kde, function(x) {thres <- threshold; ifelse(x > thres, x - thres, NA)})
# ...and reverse the raster into the polygon spatial layer
kde.delin <- rasterToPolygons(x = raster.kde, n=16, fun = function(x) { x >= threshold })
# a smaller bandwidth doesn't permit buffering.
life.circ.kde <- gUnaryUnion(kde.delin)

# ggplot figures:
g <- ggmap(sz.basemap, darken = 0.6) +
  layer_spatial(Shenzhen, alpha=0.3) +
  layer_spatial(act.cloud, color='#8DEEEE', size=0.2, alpha=0.01) +
  layer_spatial(raster.kde.cut, alpha=0.5) +
  scale_fill_scico(palette = "lapaz", na.value=NA) +
  layer_spatial(life.circ.kde %>% spTransform(CRS("+proj=longlat +datum=WGS84")), fill=NA, color = 'white') +
  annotation_north_arrow(location='tr', style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl") +
  xlab(NULL) + ylab(NULL) + guides(fill="none")

plot(g)
ggsave('./product/comparison_individual/comparison_kde_250m.png', width = 10, height=5.45)

# increase threshold, but bandwidth remain small
threshold <- 0.05
raster.kde.cut <- calc(raster.kde, function(x) {thres <- threshold; ifelse(x > thres, x - thres, NA)})
# and reverse the raster into the polygon spatial layer
kde.delin <- rasterToPolygons(x = raster.kde, n=16, fun = function(x) { x >= threshold })
# a smaller bandwidth doesn't permit buffering.
life.circ.kde <- gUnaryUnion(kde.delin)

# ggplot figures:
g <- ggmap(sz.basemap, darken = 0.6) +
  layer_spatial(Shenzhen, alpha=0.3) +
  layer_spatial(act.cloud, color='#8DEEEE', size=0.2, alpha=0.01) +
  layer_spatial(raster.kde.cut, alpha=0.5) +
  scale_fill_scico(palette = "lapaz", na.value=NA) +
  layer_spatial(life.circ.kde %>% spTransform(CRS("+proj=longlat +datum=WGS84")), fill=NA, color = 'white') +
  annotation_north_arrow(location='tr', style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl") +
  xlab(NULL) + ylab(NULL) + guides(fill="none")

plot(g)
ggsave('./product/comparison_individual/comparison_kde_250m_core.png', width = 10, height=5.45)

# SDE
ellipse <- calc_sde(points = act.cloud@coords)
ellipse.poly <- st_polygon(list(ellipse$LOCATIONS[, 2:3] %>% as.matrix()))
ellipso.poly.sf <- st_sfc(ellipse.poly, crs = "EPSG:32650")
ellipse.sf <- st_sf(geometry=ellipso.poly.sf)
# ggplot figures:
g <- ggmap(sz.basemap, darken = 0.6) +
  layer_spatial(Shenzhen, alpha=0.3) +
  layer_spatial(act.cloud, color='#8DEEEE', size=0.2, alpha=0.01) +
  layer_spatial(ellipse.sf, color='white', fill=lapaz.blue, alpha=0.5) +
  annotation_north_arrow(location='tr', style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl") +
  xlab(NULL) + ylab(NULL) + guides(fill="none")

plot(g)
ggsave('./product/comparison_individual/comparison_ellipse.png', width = 10, height=5.45)

# MCP
ch <- chull(act.cloud@coords)
cv <- act.cloud@coords[c(ch, ch[1]), ]
cvxhull <- SpatialPolygons(list(Polygons(list(Polygon(cv)), ID=1)), proj4string = crs.utm)
g <- ggmap(sz.basemap, darken = 0.6) +
  layer_spatial(Shenzhen, alpha=0.3) +
  layer_spatial(act.cloud, color='#8DEEEE', size=0.2, alpha=0.01) +
  layer_spatial(cvxhull, color='white', fill=lapaz.blue, alpha=0.5) +
  annotation_north_arrow(location='tr', style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl") +
  xlab(NULL) + ylab(NULL) + guides(fill="none")

plot(g)
ggsave('./product/comparison_individual/comparison_convexhull.png', width = 10, height=5.45)


# TDA
toporaster <- IndividAnchor[startsWith(names(IndividAnchor), '082455786')]
den.delin.list <- lapply(toporaster, rasterToPolygons, n=16)
den.delin <- Reduce(rbind, den.delin.list)
life.circ.den <- gBuffer(gUnaryUnion(den.delin), width=100)

g <- ggmap(sz.basemap, darken = 0.6) +
  layer_spatial(Shenzhen, alpha=0.3) +
  layer_spatial(act.cloud, color='#8DEEEE', size=0.2, alpha=0.01)
for(i in seq_along(toporaster)) {
  g <- g + new_scale_fill() +
    layer_spatial(toporaster[[i]], alpha=0.5) +
    scale_fill_scico(palette = "lapaz", na.value=NA, guide='none')
  }
g <- g + layer_spatial(life.circ.den, color='white', fill=NA) +
  annotation_north_arrow(location='tr', style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl") +
  xlab(NULL) + ylab(NULL) + guides(fill="none")

plot(g)
ggsave('./product/comparison_individual/comparison_topo.png', width = 10, height=5.45)



# Activity space statisics calculation
load('./image_result/new_result/fineCommuneRes_05.RData')
load('./image_result/new_result/CommuneCirc_05.RData')

# AS area in km^2
areaCalc <- function(x) {
  value.r <- reclassify(x, matrix(c(-Inf, Inf, 1, NA, NA, 0), byrow = T, nrow = 2))
  cell.cnt <- cellStats(value.r, sum)
  area <- prod(res(x) / 1e3) * cell.cnt
  return(area)
}

# descriptive statistics
sapply(CommuneCirc, areaCalc) %>% psych::describe()
sqrt(sapply(CommuneCirc, areaCalc) / pi) %>% describe()

# number of event records
sapply(location.Points.vIs, nrow) %>% psych::describe()
sapply(location.Points.Vis, nrow) %>% psych::describe()

# number of anchors
sapply(fineAnchorRes, function(x) {sum(sapply(x, function(y) length(y$repreCompo$cell)))}) %>% 
  psych::describe()

# Sketch Area -------------------------------------------------------------
