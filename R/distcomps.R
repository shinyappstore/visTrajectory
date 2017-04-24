#' Compute pairwise correlation-based distance
#' 
#' Pearson correlation distance for continuous
#' data: d = (1- cor(x_i, x_j))/2 where cor(x_i, x_j)
#' is the correlation between vector x_i and vector x_j.
#'
#' @param X A data matrix, e.g. gene expression
#' @param method a character string indicating which correlation coefficient 
#' is to be computed. One of "pearson" (default), "kendall", or 
#' "spearman": can be abbreviated.
#' @param scale A boolean indicating whether to normalize the
#' columns (samples) of the data to the even sum.
#' @param base A numeric value for the shared column sum, if scale is TRUE.  
#' @param log_trans A boolean indicating whether to log transform
#' the data prio to distance computation (log(X + 1)). Default is FALSE.
#' @param log_base A number indicating base for log transformation. 
#' Default is 10.
#' 
#' @return A dissimilarity matrix, D.
#' @export
cor_dist <- function(X, method = "pearson",
                     scale = TRUE, base = 1e6,
                     log_trans = TRUE, log_base = 10){
  if(scale) { # scale samples to the same base level
    X <- apply(X, 2, function(x) x/sum(x)*base)
  }
  if (log_trans) { # log-transform data
    logexp <- log(X + 1, base = log_base) 
  }
  # Compute correlation based distance
  corDist <- (1 - cor(X, method = method))/2
  D <- as.matrix(corDist) 
  return(D)
}

#' Compute Jaccard distance matrix
#' 
#' Jaccard (or binary) distance, implemented for
#' microbial community composition data. Include
#' filtering options before distance computation.
#'
#' @param X A data matrix, e.g. microbial counnts 
#' @param min_row_prevalence An integer indicating the minimum
#' prevalence (non-zero occurance) of a feature (species) across
#' samples to be left for distance computation.
#' @param min_row_sum A boolean indicating whether to normalize the
#' columns (samples) of the data to the even sum.
#' 
#' @return A dissimilarity matrix, D.
#' @export
jacc_dist <- function(X, min_row_sum = 100, min_row_prevalence = 5){
  X <- X[rowSums(X > 0) >= min_row_prevalence, ]
  X <- X[rowSums(X) >= min_row_sum, ]
  jaccDist <- dist(t(X), method = "binary")
  D <- as.matrix(jaccDist) 
  return(D)
}


#' Rank based transform distances to triangular distribution
#' 
#' Transforms the dissimularities to follow a triangular 
#' distribution on [0, 1] interval to better match distances
#' on scalar coordinates that are approximately uniformly
#' distributed.
#'
#' @param D A dissimilarity matrix.
#' 
#' @return A transformed dissimilarity matrix, D.
#' @export
transform_dist <- function(D) {
  dvec0 <- D[lower.tri(D)]
  dvec <- rank(dvec0)
  dvec <- dvec/max(dvec)
  dvec <- 1 - sqrt(1 - dvec)
  D <- matrix(0, nrow = nrow(D), ncol = ncol(D))
  D[lower.tri(D)] <- dvec
  D <- D + t(D)
  return(D)
}


#' k-nearest neighbours distance.
#' 
#' Compute average distance to k-nearest neighbours
#' for each data point from a dissimilarity matrix.
#' This approximates the data densities around 
#' each data point.
#' 
#' @param D A dissimilarity matrix.
#' @param K An integer indicating the number of points 
#' to include in the neighbourhood. If not specified 
#' one tenth of the total number of data points is used.
#' 
#' @return A list with mean distances to k-nearest 
#' neighbours and the associated variance of these
#' k-distances for each data point.
#' 
#' @export
kNN_dist <- function(D, K = NULL) {
  if (is.null(K) || is.na(K)) K <- floor(ncol(D)/10)
  D_kSort <- apply(D, 2, function(x) {sort(x)})
  mean <- apply(D_kSort, 2, function(x){
    mean(x[2:(K+1)]) 
  })
  var <- apply(D_kSort, 2, function(x){
    var(x[2:(K+1)]) 
  })
  return(list(mean = mean, var = var))
}


#' Reshape distance matrix to long format
#' 
#' Internal function for getting distances in a long data.frame format,
#' a format used in fit_buds function 
#'
#' @param D A square matrix of pairwise dissimilarities
#' @param ... Other parameters for princurve::principal.curve function.
#' 
#' @return A data frame to be used as input data for BUDS.
get_dist_df <- function(D, K = NULL) {
  N <- nrow(D)
  rownames(D) <- colnames(D) <- 1:N
  distDF <- reshape2::melt(D)
  colnames(distDF) <- c("j", "i", "d")
  # Keep only (n-1)n/2 distances as you assume they are symmetric
  distDF <- distDF[distDF$i < distDF$j, ]
  # Delete all the missing value distances
  distDF <- distDF[!is.na(distDF$d), ]
  # Compute average distance to kNN
  dkNN <- kNN_dist(D, K = K)
  distDF$i_sigma_K <- dkNN$mean[distDF$i]
  distDF$j_sigma_K <- dkNN$mean[distDF$j]
  distDF$max_sigma_K <- pmax(distDF$i_sigma_K, distDF$j_sigma_K)
  distDF$min_sigma_K <- pmin(distDF$i_sigma_K, distDF$j_sigma_K)
  return(distDF)
}
