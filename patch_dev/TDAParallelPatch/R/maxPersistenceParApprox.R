#' This function is similar to maxPersistencePara, yet only randomly extract a sample from the dataset
#' to estimate the confidence band. Parallel computation is also adopted.
#'
#' @param FUN the distance function adopted for filtration, such as `dtm`, `kde`, etc.
#' @param parameters a numeric vector for tuning. The optimal choice among them will be chosen.
#' @param X the dataset, a matrix.
#' @param lim the range of each dimension of the grid.
#' @param by the resolution of the grid.
#' @param maxdimension the maximum dimension to compute PH.
#' @param sublevel a logical value, whether the sublevel set PH should be computed.
#' @param library a string specifying which library should be used to compute.
#' @param B the number of bootstraps.
#' @param alpha confident level
#' @param bandFUN the function to be used in the computation, either `"bootstrapDiagram"`, `"bootstrapBand"` or `"sampleDiagram"`.
#' @param distance the same as `maxPersistence`.
#' @param dimension the same as `maxPersistence`.
#' @param p the same as `maxPersistence`.
#' @param printProgress the same as `maxPersistence`.
#' @param weight the same as `maxPersistence`.
#' @param corePermit the number of CPU cores permitted for tuning.
#' @param os the OS.
#'
#' @title Parallel Version of Maximal Persistence Method (Random Sampling Estimation)
#' @description The same function detecting maximal persistence as that of `TDA` package. In the function multiple computation kernel would be opened and parallel computation will be conducted using all of them.
#' @export

maxPersistenceParApprox=function(FUN, parameters, X, lim, by, approxSize=nrow(X), maxdimension = length(lim)/2 - 1,
                            sublevel = TRUE, library = "GUDHI", B = 30, alpha = 0.05,
                            bandFUN = "sampleBand", distance = "bottleneck", dimension = min(1, maxdimension),
                            p = 1, printProgress = FALSE, weight=NULL, corePermit = parallel::detectCores(),
                            os = 'windows')
{
  if (!is.function(FUN)) {
    stop("FUN should be a function")
  }
  if (!is.numeric(parameters)) {
    stop("parameters should be a numeric vector")
  }
  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  if (!is.numeric(lim) || length(lim)%%2 != 0) {
    stop("lim should be either a numeric matrix or a numeric vector of even elements")
  }
  if (!is.numeric(by) || any(by <= 0)) {
    stop("by should be positive")
  }
  if (2 * NCOL(X) != length(lim)) {
    stop("dimension of X does not match with lim")
  }
  if (length(by) != 1 && length(by) != NCOL(X)) {
    stop("by should be either a number or a vector of length equals dimension of grid")
  }
  if (!is.numeric(maxdimension) || length(maxdimension) != 1 || maxdimension < 0) {
    stop("maxdimnsion should be a nonnegative integer")
  }
  if (!is.logical(sublevel)) {
    stop("sublevel should be logical")
  }
  if (library == "gudhi" || library == "Gudhi") {
    library <- "GUDHI"
  }
  if (library == "dionysus" || library == "DIONYSUS") {
    library <- "Dionysus"
  }
  if (library == "phat" || library == "Phat") {
    library <- "PHAT"
  }
  if (library != "GUDHI" && library != "Dionysus" && library != "PHAT") {
    stop("library for computing persistence diagram should be a string: either 'GUDHI', 'Dionysus', or 'PHAT'")
  }
  if (!is.numeric(B) || length(B) != 1 || B < 1) {
    stop("B should be a positive integer")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha < 0 || alpha > 1) {
    stop("alpha should be a number between 0 and 1")
  }
  if (bandFUN != "sampleBand" && bandFUN != "sampleDiagram") {
    stop("bandFUN should be a string: either 'sampleBand' or 'sampleDiagram'")
  }
  if (distance != "wasserstein" && distance != "bottleneck") {
    stop("distance should be a string: either 'bottleneck' or 'wasserstein'")
  }
  if (!is.numeric(dimension) || any(dimension < 0 | dimension > maxdimension)) {
    stop("dimension should be a integer or a vector of integer, with the range between 0 and maxdimension")
  }
  if (!is.numeric(p) || length(p) != 1 || p < 1) {
    stop("p should be a positive integer")
  }
  # if (!is.logical(parallel)) {
  #   stop("parallel should be logical")
  # }
  if (!is.logical(printProgress)) {
    stop("printProgress should be logical")
  }
  if (!is.null(weight) && (!is.numeric(weight) || (length(weight) != 1 && length(weight) != NROW(X)))) {
    stop("weight should be either NULL, a number, or a vector of length equals the number of sample")
  }
  if (!is.numeric(corePermit) || length(corePermit) != 1 || corePermit < 1) {
    stop("corePermit should be a positive integer")
  }
  if (os != 'windows' && os != 'linux'){
    stop("os should be either 'windows' or 'linux'")
  }

  X <- as.matrix(X)
  maxdimension <- min(maxdimension, NCOL(X) - 1)
  Kseq <- length(parameters)
  eps <- numeric(Kseq)
  numberSignificant <- numeric(Kseq)
  significantPers <- numeric(Kseq)
  Pers <- list()
  corePermit <- min(length(parameters), corePermit)


  x.seq <- seq(lim[1,1], lim[2,1], by)
  y.seq <- seq(lim[1,2], lim[2,2], by)
  Grid <- expand.grid(x.seq, y.seq)

  maxPersKernel <- function(i){
    # calculating Pers
    if (is.null(weight)) {
      Diag <- TDA::gridDiag(X = X, FUN = FUN, lim = lim, by = by,
                            FUNvalues = NULL, maxdimension = maxdimension,
                            sublevel = sublevel, library = library, location = FALSE,
                            printProgress = FALSE, diagLimit = NULL, parameters[i])[["diagram"]]
    }
    else {
      Diag <- TDA::gridDiag(X = X, FUN = FUN, lim = lim, by = by,
                            FUNvalues = NULL, maxdimension = maxdimension,
                            sublevel = sublevel, library = library, location = FALSE,
                            printProgress = FALSE, diagLimit = NULL, weight = weight,
                            parameters[i])[["diagram"]]
    }
    Diag[1, 3] <- Diag[1, 2]
    Pers <- cbind(Diag[, 1], Diag[, 3] - Diag[, 2])
    colnames(Pers) <- c("dimension", "Persistence")
    # calculating eps, numberSignificant & significantPers
    if (bandFUN == "sampleDiagram") {
      eps <- sampleDiagram(X = X, FUN = FUN, lim = lim, by = by, sampleSize = approxSize, maxdimension = maxdimension,
                           sublevel = sublevel, library = library, B = B, alpha = alpha,
                           distance = distance, dimension = dimension, p = p, parallel = TRUE,
                           printProgress = FALSE, weight = weight, parameters[i])
      selctDim <- rep(FALSE, NROW(Pers))
      for (dim in dimension) {
        selctDim <- selctDim | (Pers[, 1] == dim)
      }
      numberSignificant <- sum(Pers[selctDim, 2] > (2 * eps))
      significantPers <- sum(pmax(0, Pers[selctDim, 2] - (2 * eps)))
    }
    else {
      eps <- sampleBand(X = X, FUN = FUN, Grid = Grid, sampleSize=approxSize,
                        B = B, alpha = alpha, parallel = TRUE, printProgress = TRUE,
                        weight = weight, parameters[i])[["width"]]
      numberSignificant <- sum(Pers[, 2] > (2 * eps))
      significantPers <- sum(pmax(0, Pers[, 2] - (2 * eps)))
    }
    if (printProgress) {
      cat("Calculation of parameter ", parameters[i], " completed.\n")
    }
    return(list(eps=eps, numberSignificant=numberSignificant, significantPers=significantPers, Pers=Pers))
  }

  if(os == 'windows')
  {
    clstrCPU <- parallel::makeCluster(corePermit)
    varlist <- c("X", "lim", "by", "maxdimension", "sublevel", "library", "parameters", "B",
                 "alpha", "distance", "dimension", "weight", "p", "Grid", "eps", "approxSize",
                 "bandFUN", "numberSignificant", "significantPers", "Pers", "printProgress", "FUN",
                 "sampleDiagram", "sampleBand")
    tryCatch({
      parallel::clusterEvalQ(clstrCPU, {library(TDA)})
      parallel::clusterExport(clstrCPU, varlist, envir = environment())
      parRes <- parallel::parLapply(clstrCPU, seq(along=parameters), maxPersKernel)
    }, error=function(e){
      stop(conditionMessage(e))
    }, finally = {parallel::stopCluster(clstrCPU)})
  }
  else{
    parRes <- parallel::mclapply(seq(along=parameters), maxPersKernel, mc.cores = corePermit)
  }

  if (printProgress) {
    cat("\n")
  }

  for(i in seq(Kseq)){
    numberSignificant[i] <- parRes[[i]]$numberSignificant
    significantPers[i] <- parRes[[i]]$significantPers
    eps[i] <- parRes[[i]]$eps
    Pers[[i]] <- parRes[[i]]$Pers
  }

  out <- list(parameters = parameters, sigNumber = numberSignificant,
              sigPersistence = significantPers, bands = eps, Persistence = Pers,
              bandFUN = bandFUN, dimension = dimension)
  class(out) <- "maxPersistence"
  attributes(out)[["call"]] <- match.call()
  return(out)
}
