# Helper function to create a lower triangular distance matrix (column-wise storage)
# from a full distance matrix
make_rdist <- function(d) {
  n <- nrow(d)
  rdist <- numeric(n * (n - 1) / 2)
  idx <- 1
  for (j in 1:(n - 1)) {
    for (i in (j + 1):n) {
      rdist[idx] <- d[i, j]
      idx <- idx + 1
    }
  }
  rdist
}

# Create test data: simple 2D points that form clear clusters
set.seed(42)
n_points <- 30
cluster1 <- matrix(rnorm(20, mean = 0, sd = 0.5), ncol = 2)
cluster2 <- matrix(rnorm(20, mean = 5, sd = 0.5), ncol = 2)
cluster3 <- matrix(rnorm(20, mean = 10, sd = 0.5), ncol = 2)
test_data <- rbind(cluster1, cluster2, cluster3)
n <- nrow(test_data)
full_dist <- as.matrix(dist(test_data))
rdist <- make_rdist(full_dist)

test_that("pam returns valid KmedoidsResult", {
  result <- pam(rdist, n, k = 3, maxiter = 100)
  
  expect_s4_class(result, "KmedoidsResult")
  expect_true(result@cost >= 0)
  expect_length(result@medoids, 3)
  expect_length(result@assignment, n)
  expect_true(all(result@medoids >= 0 & result@medoids < n))
  expect_true(all(result@assignment >= 1 & result@assignment <= 3))
})

test_that("fastpam returns valid KmedoidsResult", {
  result <- fastpam(rdist, n, k = 3, maxiter = 100)
  
  expect_s4_class(result, "KmedoidsResult")
  expect_true(result@cost >= 0)
  expect_length(result@medoids, 3)
  expect_length(result@assignment, n)
  expect_true(all(result@medoids >= 0 & result@medoids < n))
  expect_true(all(result@assignment >= 1 & result@assignment <= 3))
})

test_that("fastpam works with BUILD initializer", {
  result <- fastpam(rdist, n, k = 3, maxiter = 100, initializer = "BUILD")
  
  expect_s4_class(result, "KmedoidsResult")
  expect_true(result@cost >= 0)
  expect_length(result@medoids, 3)
})

test_that("fastpam works with LAB initializer", {
  result <- fastpam(rdist, n, k = 3, maxiter = 100, initializer = "LAB")
  
  expect_s4_class(result, "KmedoidsResult")
  expect_true(result@cost >= 0)
  expect_length(result@medoids, 3)
})

test_that("fastclara returns valid KmedoidsResult", {
  result <- fastclara(rdist, n, k = 3, maxiter = 100)
  
  expect_s4_class(result, "KmedoidsResult")
  expect_true(result@cost >= 0)
  expect_length(result@medoids, 3)
  expect_length(result@assignment, n)
  expect_true(all(result@medoids >= 0 & result@medoids < n))
  expect_true(all(result@assignment >= 1 & result@assignment <= 3))
})

test_that("fastclara works with different sampling parameters", {
  result <- fastclara(rdist, n, k = 3, numsamples = 3, sampling = 0.5)
  
  expect_s4_class(result, "KmedoidsResult")
  expect_true(result@cost >= 0)
  expect_length(result@medoids, 3)
})

test_that("fastclarans returns valid KmedoidsResult", {
  result <- fastclarans(rdist, n, k = 3)
  
  expect_s4_class(result, "KmedoidsResult")
  expect_true(result@cost >= 0)
  expect_length(result@medoids, 3)
  expect_length(result@assignment, n)
  expect_true(all(result@medoids >= 0 & result@medoids < n))
  expect_true(all(result@assignment >= 1 & result@assignment <= 3))
})

test_that("fastclarans works with different parameters", {
  result <- fastclarans(rdist, n, k = 3, numlocal = 3, maxneighbor = 0.05)
  
  expect_s4_class(result, "KmedoidsResult")
  expect_true(result@cost >= 0)
  expect_length(result@medoids, 3)
})

test_that("all algorithms produce unique medoids", {
  result_pam <- pam(rdist, n, k = 3, maxiter = 100)
  result_fastpam <- fastpam(rdist, n, k = 3, maxiter = 100)
  result_fastclara <- fastclara(rdist, n, k = 3, maxiter = 100)
  result_fastclarans <- fastclarans(rdist, n, k = 3)
  

  expect_equal(length(unique(result_pam@medoids)), 3)
  expect_equal(length(unique(result_fastpam@medoids)), 3)
  expect_equal(length(unique(result_fastclara@medoids)), 3)
  expect_equal(length(unique(result_fastclarans@medoids)), 3)
})

test_that("seed parameter produces reproducible results", {
  result1 <- fastpam(rdist, n, k = 3, seed = 42)
  result2 <- fastpam(rdist, n, k = 3, seed = 42)
  
  expect_equal(result1@medoids, result2@medoids)
  expect_equal(result1@assignment, result2@assignment)
  expect_equal(result1@cost, result2@cost)
})

test_that("different seeds can produce different results",
{
  result1 <- fastclarans(rdist, n, k = 3, seed = 1)
  result2 <- fastclarans(rdist, n, k = 3, seed = 999)
  
  # Results may or may not differ, but both should be valid
  expect_s4_class(result1, "KmedoidsResult")
  expect_s4_class(result2, "KmedoidsResult")
})

test_that("k=1 produces single cluster", {
  result <- fastpam(rdist, n, k = 1, maxiter = 100)
  
  expect_length(result@medoids, 1)
  expect_true(all(result@assignment == 1))
})

test_that("k=2 produces two clusters", {
  result <- fastpam(rdist, n, k = 2, maxiter = 100)
  
  expect_length(result@medoids, 2)
  expect_true(all(result@assignment %in% c(1, 2)))
})

# Test with smaller dataset
test_that("algorithms work with small dataset", {
  small_data <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  small_n <- nrow(small_data)
  small_dist <- as.matrix(dist(small_data))
  small_rdist <- make_rdist(small_dist)
  

  result_pam <- pam(small_rdist, small_n, k = 2, maxiter = 100)
  result_fastpam <- fastpam(small_rdist, small_n, k = 2, maxiter = 100)
  
  expect_s4_class(result_pam, "KmedoidsResult")
  expect_s4_class(result_fastpam, "KmedoidsResult")
  expect_length(result_pam@medoids, 2)
  expect_length(result_fastpam@medoids, 2)
})
