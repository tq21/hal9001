test_that("cluster_knot_screening using pam/kmeans works", {
  # d = 1, continuous
  X <- matrix(rnorm(100), ncol = 1)
  kmeans_basis_list <- cluster_knot_screening(X, "kmeans", 10)
  pam_basis_list <- cluster_knot_screening(X, "pam", 10)
  expect_equal(length(kmeans_basis_list), 10)
  expect_equal(length(pam_basis_list), 10)

  # d = 1, binary
  X <- matrix(rep(c(1, 0), 50), ncol = 1)
  kmeans_basis_list <- cluster_knot_screening(X, "kmeans", 10)
  pam_basis_list <- cluster_knot_screening(X, "pam", 10)
  expect_equal(length(kmeans_basis_list), 2)
  expect_equal(length(pam_basis_list), 2)
})

test_that("cluster_knot_screening returns the correct number of knots
          when num_clusters is unspecified", {
  X <- matrix(rnorm(300), ncol = 3)
  default_num_clusters <- nrow(X) / 2
  num_combos <- choose(3, 1) + choose(3, 2) + choose(3, 3)
  kmeans_basis_list <- cluster_knot_screening(X, "kmeans")
  pam_basis_list <- cluster_knot_screening(X, "pam")
  expect_equal(length(kmeans_basis_list), default_num_clusters*num_combos)
  expect_equal(length(pam_basis_list), default_num_clusters*num_combos)
})

test_that("cluster_knot_screening returns the correct number of knots
          when num_clusters is specified", {
  X <- matrix(rnorm(300), ncol = 3)
  num_clusters <- 10
  num_combos <- choose(3, 1) + choose(3, 2) + choose(3, 3)
  kmeans_basis_list <- cluster_knot_screening(X, "kmeans", num_clusters)
  pam_basis_list <- cluster_knot_screening(X, "pam", num_clusters)
  expect_equal(length(kmeans_basis_list), num_clusters*num_combos)
  expect_equal(length(pam_basis_list), num_clusters*num_combos)
})

test_that("cluster_knot_screening returns the correct number of knots
          when max_degree is specified", {
  X <- matrix(rnorm(300), ncol = 3)
  max_degree <- 2
  num_combos <- choose(3, 1) + choose(3, 2)
  kmeans_basis_list <- cluster_knot_screening(X, "kmeans", max_degree = max_degree)
  pam_basis_list <- cluster_knot_screening(X, "pam", max_degree = max_degree)
  expect_equal(length(kmeans_basis_list), nrow(X)/2*num_combos)
  expect_equal(length(pam_basis_list), nrow(X)/2*num_combos)
})

test_that("cluster_knot_screening returns the correct number of knots when the
          number of unique values in X is less than the number of clusters", {
  X <- matrix(rep(c(1, 2, 3), 50), ncol = 1)
  kmeans_basis_list <- cluster_knot_screening(X, "kmeans", 10)
  pam_basis_list <- cluster_knot_screening(X, "pam", 10)
  expect_equal(length(kmeans_basis_list), 3)
  expect_equal(length(pam_basis_list), 3)
})

test_that("cluster_knot_screening gives error when algorithm is neither
          pam nor kmeans", {
  X <- matrix(rnorm(300), ncol = 3)
  expect_error(cluster_knot_screening(X, "random"))
})
