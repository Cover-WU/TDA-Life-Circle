#' This function is the confidence estimator similar to the `TDA::bootstrapDiagram`,
#' usually for big dataset where full size analysis is computationally intractable.
#'
#' @param X the dataset, a matrix.
#' @param FUN the distance function adopted for filtration, such as `dtm`, `kde`, etc.
#' @param lim the range of each dimension of the grid.
#' @param by the resolution of the grid.
#' @param sampleSize the sample size of the X during sampling step.
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
#' @param weight the same as  `maxPersistence`.
#'
#' @title Sampling Confidence Set for a Persistence Diagram, using the Bottleneck Distance (or the Wasserstein distance).
#' @export

sampleDiagram = function (X, FUN, lim, by, sampleSize=nrow(X), maxdimension = length(lim)/2 - 1,
          sublevel = TRUE, library = "GUDHI", B = 30, alpha = 0.05,
          distance = "bottleneck", dimension = min(1, maxdimension),
          p = 1, parallel = FALSE, printProgress = FALSE, weight=NULL, ...)
{
  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  if (!is.numeric(lim) || length(lim)%%2 != 0) {
    stop("lim should be either a numeric matrix or a numeric vector of even elements")
  }
  if (!is.numeric(by) || any(by <= 0)) {
    stop("by should be positive")
  }
  if (!is.numeric(sampleSize) || sampleSize > NROW(X) || sampleSize < 0) {
    stop("sampleSize should be a positive number not greater than row number of X")
  }
  if (2 * NCOL(X) != length(lim)) {
    stop("dimension of X does not match with lim")
  }
  if (length(by) != 1 && length(by) != NCOL(X)) {
    stop("by should be either a number or a vector of length equals dimension of grid")
  }
  if (!is.numeric(maxdimension) || length(maxdimension) !=
      1 || maxdimension < 0) {
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
  if (library != "GUDHI" && library != "Dionysus" && library !=
      "PHAT") {
    stop("library for computing persistence diagram should be a string: either 'GUDHI', 'Dionysus', or 'PHAT'")
  }
  if (!is.numeric(B) || length(B) != 1 || B < 1) {
    stop("B should be a positive integer")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <
      0 || alpha > 1) {
    stop("alpha should be a number between 0 and 1")
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
  if (!is.logical(parallel)) {
    stop("parallel should be logical")
  }
  if (!is.logical(printProgress)) {
    stop("printProgress should be logical")
  }
  if (!is.null(weight) && (!is.numeric(weight) || (length(weight) !=
                                                   1 && length(weight) != NROW(X)))) {
    stop("weight should be either NULL, a number, or a vector of length equals the number of sample")
  }
  X <- as.matrix(X)
  maxdimension <- min(maxdimension, NCOL(X) - 1)

  if (is.null(weight)){
    Diag <- TDA::gridDiag(X = X, FUN = FUN, lim = lim, by = by,
                       FUNvalues = NULL, maxdimension = maxdimension, sublevel = sublevel,
                       library = library, location = FALSE, printProgress = FALSE,
                       diagLimit = NULL, ...)[["diagram"]]
    if (distance == "wasserstein") {
      boostFUN <- function(i) {
        I <- sample(NROW(X), replace = FALSE, size = sampleSize)
        Diag1 <- TDA::gridDiag(X = X[I, , drop = FALSE],
                          FUN = FUN, lim = lim, by = by, FUNvalues = NULL,
                          maxdimension = maxdimension, sublevel = sublevel,
                          library = library, location = FALSE, printProgress = FALSE,
                          diagLimit = NULL, ...)[["diagram"]]
        width1 <- TDA::wasserstein(Diag, Diag1, p = p, dimension = dimension)
        if (printProgress) {
          cat(i, " ")
        }
        return(width1)
      }
    }
    else {
      boostFUN <- function(i) {
        I <- sample(NROW(X), replace = FALSE, size = sampleSize)
        Diag1 <- TDA::gridDiag(X = X[I, , drop = FALSE],
                          FUN = FUN, lim = lim, by = by, FUNvalues = NULL,
                          maxdimension = maxdimension, sublevel = sublevel,
                          library = library, location = FALSE, printProgress = FALSE,
                          diagLimit = NULL, ...)[["diagram"]]
        width1 <- TDA::bottleneck(Diag, Diag1, dimension = dimension)
        if (printProgress) {
          cat(i, " ")
        }
        return(width1)
      }
    }
  }
  else {
    Diag <- gridDiag(X = X, FUN = FUN, lim = lim, by = by,
                     FUNvalues = NULL, maxdimension = maxdimension, sublevel = sublevel,
                     library = library, location = FALSE, printProgress = FALSE,
                     diagLimit = NULL, weight = weight, ...)[["diagram"]]
    if (distance == "wasserstein") {
      boostFUN <- function(i) {
        size = sum(weight); prob = weight
        weightBoost <- rowSums(stats::rmultinom(n = size %/% (2 ^ 31 - 1), size = 2 ^ 31 - 1, prob = prob)) +
          as.vector(stats::rmultinom(n = 1, size = size %% (2^31-1), prob = prob))

        I <- sample(NROW(X), replace = FALSE, size = sampleSize)
        weightBoost <- ifelse(length(weightBoost) > 1, weightBoost[I], weightBoost)
        Diag1 <- gridDiag(X = X[I, , drop = FALSE], FUN = FUN, lim = lim,
                          by = by, FUNvalues = NULL, maxdimension = maxdimension,
                          sublevel = sublevel, library = library, location = FALSE,
                          printProgress = FALSE, diagLimit = NULL, weight = weightBoost,
                          ...)[["diagram"]]
        width1 <- wasserstein(Diag, Diag1, p = p, dimension = dimension)
        if (printProgress) {
          cat(i, " ")
        }
        return(width1)
      }
    }
    else {
      boostFUN <- function(i) {
        size = sum(weight); prob = weight
        weightBoost <- rowSums(stats::rmultinom(n = size %/% (2 ^ 31 - 1), size = 2 ^ 31 - 1, prob = prob)) +
          as.vector(stats::rmultinom(n = 1, size = size %% (2^31-1), prob = prob))

        I <- sample(NROW(X), replace = FALSE, size = sampleSize)
        weightBoost <- ifelse(length(weightBoost) > 1, weightBoost[I], weightBoost)
        Diag1 <- gridDiag(X = X[I, , drop = FALSE], FUN = FUN, lim = lim,
                          by = by, FUNvalues = NULL, maxdimension = maxdimension,
                          sublevel = sublevel, library = library, location = FALSE,
                          printProgress = FALSE, diagLimit = NULL, weight = weightBoost,
                          ...)[["diagram"]]
        width1 <- bottleneck(Diag, Diag1, dimension = dimension)
        if (printProgress) {
          cat(i, " ")
        }
        return(width1)
      }
    }
  }

  if (parallel) {
    boostLapply <- parallel::mclapply
  }
  else {
    boostLapply <- lapply
  }
  if (printProgress) {
    cat("Bootstrap: ")
  }
  width <- boostLapply(seq_len(B), FUN = boostFUN)
  if (printProgress) {
    cat("\n")
  }
  width <- stats::quantile(unlist(width), 1 - alpha)
  return(width)
}
