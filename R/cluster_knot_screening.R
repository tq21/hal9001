#' Function to screen down knot-points using clustering algorithms
#'
#' @description Screen down HAL knot-points using clustering algorithms to
#' reduce the subsequent computational and memory cost of making the augmented
#' design matrix and fitting HAL.
#'
#' @param X An input `matrix` containing observations and covariates following
#' standard conventions in problems of statistical learning.
#' @param algorithm The clustering algorithm to be used. Currently, only
#' \code{pam} and \code{kmeans} are supported. See \section{Details} for more
#' information.
#' @param num_clusters The number of desired clusters (knot-points) for each
#' HAL basis function.
#' @param max_degree The highest order of interaction terms for which the basis
#' functions ought to be generated. The default (NULL) corresponds to generating
#' basis functions for the full dimensionality of the input matrix.
#'
#' @importFrom stats kmeans
#' @importFrom cluster pam
#' @importFrom data.table is.data.table as.data.table uniqueN
#'
#' @returns A `list` of basis functions generated for all covariates and
#' interaction thereof up to a pre-specified degree.
#'
#' @export
#'
#' @examples
#' gendata <- function(n) {
#'   W1 <- runif(n, -3, 3)
#'   W2 <- rnorm(n)
#'   W3 <- runif(n)
#'   W4 <- rnorm(n)
#'   g0 <- plogis(0.5 * (-0.8 * W1 + 0.39 * W2 + 0.08 * W3 - 0.12 * W4))
#'   A <- rbinom(n, 1, g0)
#'   Q0 <- plogis(0.15 * (2 * A + 2 * A * W1 + 6 * A * W3 * W4 - 3))
#'   Y <- rbinom(n, 1, Q0)
#'   data.frame(A, W1, W2, W3, W4, Y)
#' }
#' set.seed(1234)
#' data <- gendata(100)
#' covars <- setdiff(names(data), "Y")
#' X <- as.matrix(data[, covars, drop = FALSE])
#' screened_basis_list <- cluster_knot_screening(X, "kmeans", 10)
cluster_knot_screening <- function(X,
                                   algorithm,
                                   num_clusters = NULL,
                                   max_degree = NULL) {

  if (is.null(num_clusters)) num_clusters <- ceiling(nrow(X)^(1/3)*5)
  if (is.null(max_degree)) max_degree <- ncol(X)

  if (!is.data.table(X)) {
    X <- as.data.table(X)
  }

  # generate column index sets up to max_degree interactions
  var_idx_num <- unlist(lapply(1:min(3, max_degree), function(g) {
    combn(1:ncol(X), g, simplify = FALSE)
  }), recursive = FALSE)

  # make basis list for each column index set
  # screen down the number of knot-points using clustering algorithms
  basis_list <- unlist(lapply(var_idx_num, function(var_idx) {
    X_sub <- X[, ..var_idx]

    if (algorithm == "pam") {
      pam_obj <- pam(x = X_sub, k = min(num_clusters, uniqueN(X_sub)),
                     metric = "euclidean")
      centroids <- pam_obj$medoids
    } else if (algorithm == "kmeans") {
      k_means <- kmeans(X_sub, centers = min(num_clusters, uniqueN(X_sub)),
                        algorithm = "MacQueen")
      centroids <- k_means$centers
    } else {
      stop("knot-point screening algorithm must be either pam or kmeans")
    }

    basis <- unname(apply(centroids, 1, function(centroid) {
      return(list(cols = var_idx,
                  cutoffs = as.numeric(centroid),
                  orders = rep(0, length(var_idx))))
    }))

    return(basis)
  }), recursive = FALSE)

  return(basis_list)
}
