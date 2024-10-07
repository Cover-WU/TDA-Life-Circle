# sampleDiagram(X, dtm, lim, by, sampleSize=100, maxdimension = 2, B=B, alpha=alpha, dimension=0, m0=0.03)
# bootstrapDiagram(X, dtm, lim, by, maxdimension = 1, B=B, alpha = alpha, dimension = 0, m0=.03)


tuneParamM0test <- function(X, lim, by, grid, n.split, B, alpha, FUN=dtm, minimum=NA, maximum=0.8,
                              fluctLimit=0.1, approx=FALSE, approxSize=1e3){
  # initialization
  fluct <- Inf; smplsz.diff <- Inf; eps <- 1e-6
  start <- 0; end <- maximum; width <- end - start
  persisLen <- 0

  n <- nrow(X)
  if(is.na(minimum)){
    minimum <- 1 / n
  }

  while ((fluct >= fluctLimit * persisLen - eps || width >= 0.05 - eps) && smplsz.diff > 0) {
    m0.v <- seq(start, end, length.out=n.split)
    # print the logs
    cat("Tuning the parameters below: \n", m0.v, '\n')
    # tuning the parameter
    if(n > approxSize && approx){
      maxPers <- maxPersistenceParApprox(FUN=FUN, ifelse(m0.v < minimum, minimum, m0.v), X, lim, by, approxSize, bandFUN = "sampleDiagram", B=B, alpha = alpha, distance = "bottleneck", dimension = c(0,1), printProgress = TRUE, corePermit=n.split)
    }
    else {
      maxPers <- maxPersistencePara(FUN=FUN, ifelse(m0.v < minimum, minimum, m0.v), X, lim, by, Grid = grid, bandFUN = "bootstrapDiagram", B=B, alpha = alpha, distance = "bottleneck", dimension = c(0,1), printProgress = TRUE, corePermit=n.split)
    }
    # select the optima
    train.frame <- tibble(maxPers$parameters, maxPers$sigNumber, maxPers$sigPersistence) %>%
      set_names(c("m0", "sigNumber", "sigPersistence"))
    m.o <- train.frame %>% arrange(desc(sigNumber), desc(sigPersistence)) %>% slice(1) %>% .$m0
    maxIdx <- which(ifelse(m0.v == 0, minimum, m0.v) == m.o)
    # update the gauge
    persisLen <- maxPers$sigPersistence[maxIdx]
    cutIdx <- c(ifelse(maxIdx - 1 > 0, maxIdx - 1, 1), ifelse(maxIdx + 1 > n.split, n.split, maxIdx + 1))
    # cut the interval again for the next iteration
    start <- m0.v[cutIdx[1]]; end <- m0.v[cutIdx[2]]; width <- end - start
    # verifying the current selection
    fluct <- persisLen - mean(maxPers$sigPersistence[cutIdx])
    # compute the difference of the sample size
    start.adj <- ifelse(start == 0, minimum, start)
    smplsz.diff <- ceiling(end * n) - ceiling(start.adj * n)
  }

  result <- list(completeRes = maxPers, simpleRes = train.frame, optimal = m.o)
  return(result)
}


# ans <- tuneParamM0test(X, lim, by, grid, n.split = nNodes, B = B, alpha = alpha)

ans <- tuneParamM0test(X, lim, by, grid, n.split = nNodes, B = B, alpha = alpha, approx = T, approxSize = 1e2)
