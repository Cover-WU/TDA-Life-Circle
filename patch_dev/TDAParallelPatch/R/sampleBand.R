#' This function is the confidence estimator similar to the `TDA::bootstrapBand`,
#' usually for big dataset where full size analysis is computationally intractable.
#'
#' @param X the dataset, a matrix.
#' @param FUN the distance function adopted for filtration, such as `dtm`, `kde`, etc.
#' @param Grid the meshgrid used for computation.
#' @param sampleSize the sample size of the X during sampling step.
#' @param B the number of bootstraps.
#' @param alpha confident level
#' @param printProgress the same as `maxPersistence`.
#' @param weight the same as  `maxPersistence`.
#'
#' @title Sampling Confidence Band
#' @export


sampleBand = function (X, FUN, Grid, sampleSize, B = 30, alpha = 0.05, parallel = FALSE,
          printProgress = FALSE, weight = NULL, ...)
{
  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  if (!is.numeric(Grid) && !is.data.frame(Grid)) {
    stop("Grid should be a matrix of coordinates")
  }
  if (!is.numeric(sampleSize) || sampleSize > NROW(X) || sampleSize < 0) {
    stop("sampleSize should be a positive number not greater than row number of X")
  }
  if (!is.numeric(B) || length(B) != 1 || B < 1) {
    stop("B should be a positive integer")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <
      0 || alpha > 1) {
    stop("alpha should be a number between 0 and 1")
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
  if (is.null(weight)) {
    ff <- FUN(X, Grid, ...)
    boostFUN <- function(i) {
      I <- sample(NROW(X), replace = FALSE, size = sampleSize)
      bootF <- FUN(X[I, , drop = FALSE], Grid, ...)
      width1 <- max(abs(ff - bootF))
      if (printProgress) {
        cat(i, " ")
      }
      return(width1)
    }
  }
  else {
    ff <- FUN(X, Grid, weight = weight, ...)
    boostFUN <- function(i) {
      I <- sample(NROW(X), replace = FALSE, size = sampleSize)
      size = sum(weight); prob = weight
      weightBoost <- rowSums(stats::rmultinom(n = size %/% (2 ^ 31 - 1), size = 2 ^ 31 - 1, prob = prob)) +
        as.vector(stats::rmultinom(n = 1, size = size %% (2^31-1), prob = prob))
      weightBoost <- ifelse(length(weightBoost) > 1, weightBoost[I], weightBoost)

      bootF <- FUN(X[I, , drop = FALSE], Grid, weight = weightBoost,
                   ...)
      width1 <- max(abs(ff - bootF))
      if (printProgress) {
        cat(i, " ")
      }
      return(width1)
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
  UPband <- ff + width
  LOWband <- ff - width
  Band <- cbind(LOWband, UPband)
  return(list(width = width, fun = ff, band = Band))
}
