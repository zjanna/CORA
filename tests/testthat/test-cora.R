test_that("CORA reduces to BH under independence", {
  set.seed(42)
  m <- 100
  rawp <- runif(m)
  cor_mat <- diag(m)  # identity = no correlation

  adj_bh   <- adjust_BH(rawp)
  adj_cora <- adjust_CORA(rawp, cor_mat)

  expect_equal(adj_cora, adj_bh, tolerance = 1e-10)
})

test_that("CORA is more conservative than BH under correlation", {
  set.seed(42)
  m <- 50
  rawp <- c(runif(10, 0, 0.01), runif(40, 0, 1))

  # Block correlation
  cor_mat <- diag(m)
  cor_mat[1:10, 1:10] <- 0.5
  diag(cor_mat) <- 1

  adj_bh   <- adjust_BH(rawp)
  adj_cora <- adjust_CORA(rawp, cor_mat)

  # CORA adjusted p-values >= BH adjusted p-values
  expect_true(all(adj_cora >= adj_bh - 1e-10))

  # CORA rejects fewer or equal hypotheses
  expect_lte(sum(adj_cora < 0.05), sum(adj_bh < 0.05))
})

test_that("CORA rejections are subset of BH rejections", {
  set.seed(123)
  m <- 200
  rawp <- c(runif(30, 0, 0.005), runif(170, 0, 1))

  cor_mat <- diag(m)
  cor_mat[1:30, 1:30] <- 0.6
  diag(cor_mat) <- 1

  adj_bh   <- adjust_BH(rawp)
  adj_cora <- adjust_CORA(rawp, cor_mat)

  rej_bh   <- which(adj_bh   < 0.05)
  rej_cora <- which(adj_cora < 0.05)

  # Every CORA rejection is also a BH rejection
  expect_true(all(rej_cora %in% rej_bh))
})

test_that("adjust_CORA validates inputs", {
  expect_error(adjust_CORA("abc", diag(3)), "numeric")
  expect_error(adjust_CORA(c(0.1, NA, 0.3), diag(3)), "NA")
  expect_error(adjust_CORA(c(0.1, 1.5, 0.3), diag(3)), "between 0 and 1")
  expect_error(adjust_CORA(c(0.1, 0.2), diag(3)), "m x m")
})

test_that("all methods return valid p-values", {
  set.seed(42)
  rawp <- runif(50)

  for (fn in list(adjust_BH, adjust_BY, adjust_ABH)) {
    adj <- fn(rawp)
    expect_equal(length(adj), 50)
    expect_true(all(adj >= 0 & adj <= 1))
  }

  adj_tsbh <- adjust_TSBH(rawp, 0.05)
  expect_equal(length(adj_tsbh), 50)
  expect_true(all(adj_tsbh >= 0 & adj_tsbh <= 1))
})

test_that("cora_compare returns correct structure", {
  set.seed(42)
  m <- 100
  rawp <- runif(m)
  cor_mat <- diag(m)

  result <- cora_compare(rawp, cor_mat)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 5)  # BH, BY, ABH, TSBH, CORA
  expect_true("CORA" %in% result$method)
})

test_that("cora_compare works without cor_matrix", {
  rawp <- runif(50)
  result <- cora_compare(rawp)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 4)  # BH, BY, ABH, TSBH (no CORA)
  expect_false("CORA" %in% result$method)
})

test_that("cora_summary returns all components", {
  set.seed(42)
  m <- 100
  rawp <- c(runif(20, 0, 0.005), runif(80, 0, 1))
  names(rawp) <- paste0("gene", 1:m)
  cor_mat <- diag(m)
  cor_mat[1:20, 1:20] <- 0.4
  diag(cor_mat) <- 1

  res <- cora_summary(rawp, cor_mat)

  expect_type(res, "list")
  expect_true("comparison" %in% names(res))
  expect_true("reduction"  %in% names(res))
  expect_true("top_genes"  %in% names(res))
  expect_true("overlap"    %in% names(res))
  expect_true(res$reduction >= 0)
})

test_that("stronger correlation produces stronger reduction", {
  set.seed(42)
  m <- 200
  rawp <- c(runif(40, 0, 0.005), runif(160, 0, 1))

  # Weak correlation
  cor_weak <- diag(m)
  cor_weak[1:40, 1:40] <- 0.1
  diag(cor_weak) <- 1

  # Strong correlation
  cor_strong <- diag(m)
  cor_strong[1:40, 1:40] <- 0.7
  diag(cor_strong) <- 1

  n_weak   <- sum(adjust_CORA(rawp, cor_weak)   < 0.05)
  n_strong <- sum(adjust_CORA(rawp, cor_strong) < 0.05)

  # Stronger correlation should produce fewer rejections
  expect_lte(n_strong, n_weak)
})
