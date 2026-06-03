#' CORA: COrrelation-Redundancy-Aware FDR Adjustment
#'
#' Adjusts p-values using the CORA method, which modifies the Benjamini-Hochberg
#' procedure by incorporating pairwise gene correlations. Genes that are strongly
#' correlated with already-rejected genes receive a stricter threshold, reducing
#' redundancy in the list of differentially expressed genes.
#'
#' @param rawp Numeric vector of raw (unadjusted) p-values.
#' @param cor_matrix Square matrix of absolute pairwise correlations between
#'   genes. Must have \code{length(rawp)} rows and columns. Should contain
#'   absolute values (|r|), not signed correlations.
#'
#' @return Numeric vector of CORA-adjusted p-values (same length as \code{rawp}).
#'
#' @details
#' The CORA correction coefficient for gene at rank \eqn{i} is:
#' \deqn{c_i = \sum_{k=1}^{i-1} |r_{\sigma(k), \sigma(i)}|}
#' where \eqn{\sigma} is the ranking permutation (by ascending p-value).
#'
#' The adjusted p-values are:
#' \deqn{\tilde{p}_{(i)}^{CORA} = \min_{j \geq i} \left\{\frac{m}{j - c_j} \cdot p_{(j)}\right\}}
#'
#' When \eqn{c_i = 0} (no correlation with previous rejections), CORA reduces to BH.
#' When \eqn{c_i > 0}, CORA is more conservative than BH.
#'
#' @references
#' Zyprych-Walczak J. CORA: a COrrelation-Redundancy-Aware FDR adjustment
#' for genomic data. BMC Bioinformatics. (submitted).
#'
#' Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical
#' and powerful approach to multiple testing. J R Stat Soc B. 1995;57:289-300.
#'
#' @examples
#' # Simulate correlated p-values
#' set.seed(42)
#' m <- 100
#' rawp <- c(runif(20, 0, 0.01), runif(80, 0, 1))  # 20 true DEGs
#'
#' # Create a block-correlated structure
#' cor_mat <- diag(m)
#' cor_mat[1:20, 1:20] <- 0.5  # DEGs are correlated
#' diag(cor_mat) <- 1
#'
#' # Compare BH and CORA
#' adj_bh   <- adjust_BH(rawp)
#' adj_cora <- adjust_CORA(rawp, cor_mat)
#'
#' cat("BH rejections:",   sum(adj_bh   < 0.05), "\n")
#' cat("CORA rejections:", sum(adj_cora < 0.05), "\n")
#'
#' @export
adjust_CORA <- function(rawp, cor_matrix) {

  # Input validation
  if (!is.numeric(rawp)) stop("rawp must be a numeric vector")
  if (any(is.na(rawp))) stop("rawp must not contain NA values")
  if (any(rawp < 0 | rawp > 1)) stop("rawp values must be between 0 and 1")

  m <- length(rawp)

  if (!is.matrix(cor_matrix)) stop("cor_matrix must be a matrix")
  if (nrow(cor_matrix) != m || ncol(cor_matrix) != m)
    stop("cor_matrix must be m x m where m = length(rawp)")

  # Sort p-values
  index <- order(rawp)
  spval <- rawp[index]

  # Compute adjusted p-values (backward pass with min-envelope)
  tmp <- c(spval, 1)  # sentinel value at position m+1

  for (i in m:1) {
    if (i == 1) {
      # First gene: no previous rejections, same as BH
      tmp[i] <- min(tmp[i + 1], min((m / i) * spval[i], 1, na.rm = TRUE),
                    na.rm = TRUE)
    } else {
      # Correlation coefficient: sum of |r| with all previously ranked genes
      coeff <- sum(cor_matrix[index[1:(i - 1)], index[i]])

      if ((i - coeff) > 0) {
        tmp[i] <- min(tmp[i + 1],
                      min((m / (i - coeff)) * spval[i], 1, na.rm = TRUE),
                      na.rm = TRUE)
      } else {
        # Degenerate case: coeff >= i (extreme correlation)
        tmp[i] <- 1
      }
    }
  }

  # Map back to original order
  result <- numeric(m)
  result[index] <- tmp[1:m]
  return(result)
}


#' Benjamini-Hochberg (BH) p-value adjustment
#'
#' Standard BH step-up procedure for FDR control under independence or PRDS.
#' Wrapper around \code{\link[stats]{p.adjust}} from base R.
#'
#' @param rawp Numeric vector of raw p-values.
#' @return Numeric vector of BH-adjusted p-values.
#'
#' @references
#' Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical
#' and powerful approach to multiple testing. J R Stat Soc B. 1995;57:289-300.
#'
#' @examples
#' rawp <- c(0.001, 0.01, 0.03, 0.5, 0.9)
#' adjust_BH(rawp)
#'
#' @export
adjust_BH <- function(rawp) {
  stats::p.adjust(rawp, method = "BH")
}


#' Benjamini-Yekutieli (BY) p-value adjustment
#'
#' BY procedure for FDR control under arbitrary dependence.
#' More conservative than BH due to the harmonic sum correction factor.
#' Wrapper around \code{\link[stats]{p.adjust}} from base R.
#'
#' @param rawp Numeric vector of raw p-values.
#' @return Numeric vector of BY-adjusted p-values.
#'
#' @references
#' Benjamini Y, Yekutieli D. The control of the false discovery rate in
#' multiple testing under dependency. Ann Stat. 2001;29:1165-1188.
#'
#' @export
adjust_BY <- function(rawp) {
  stats::p.adjust(rawp, method = "BY")
}


#' Adaptive Benjamini-Hochberg (ABH) p-value adjustment
#'
#' ABH procedure that estimates the number of true nulls (m0) adaptively
#' from the data and applies BH with the estimated m0.
#'
#' @param rawp Numeric vector of raw p-values.
#' @return Numeric vector of ABH-adjusted p-values.
#'
#' @details
#' The estimator of m0 follows Benjamini and Hochberg (2000): for each
#' ordered p-value \eqn{p_{(k)}}, compute \eqn{(m+1-k)/(1-p_{(k)})} and
#' find the first index where this sequence starts increasing. When no
#' such index exists (monotonically decreasing sequence), ABH defaults
#' to \eqn{\hat{m}_0 = m}, which makes ABH identical to BH. Note that
#' \code{multtest::mt.rawp2adjp} returns NA in this edge case; this
#' implementation handles it gracefully.
#'
#' @references
#' Benjamini Y, Hochberg Y. On the adaptive control of the false discovery
#' rate in multiple testing with independent statistics. J Educ Behav Stat.
#' 2000;25:60-83.
#'
#' @export
adjust_ABH <- function(rawp) {
  m <- length(rawp)
  index <- order(rawp)
  spval <- rawp[index]
  tmp <- spval

  # Estimate m0
  h0.m <- sapply(1:m, function(k) {
    if (spval[k] < 1) (m + 1 - k) / (1 - spval[k]) else m
  })
  grab <- which(diff(h0.m) > 0)
  h0.ABH <- if (length(grab) > 0) ceiling(min(h0.m[min(grab)], m)) else m

  # BH adjusted p-values scaled by m0/m
  for (i in (m - 1):1)
    tmp[i] <- min(tmp[i + 1], min((m / i) * spval[i], 1, na.rm = TRUE),
                  na.rm = TRUE)
  tmp <- pmin(tmp * h0.ABH / m, 1)

  result <- numeric(m)
  result[index] <- tmp
  return(result)
}


#' Two-Step Benjamini-Hochberg (TSBH) p-value adjustment
#'
#' TSBH procedure that estimates the number of true nulls via a two-step
#' approach: first applies BH at level alpha/(1+alpha) to count rejections,
#' then re-applies BH at the inflated level alpha*m/m0.
#'
#' @param rawp Numeric vector of raw p-values.
#' @param alpha Significance level (default 0.05).
#' @return Numeric vector of TSBH-adjusted p-values.
#'
#' @references
#' Benjamini Y, Krieger AM, Yekutieli D. Adaptive linear step-up procedures
#' that control the false discovery rate. Biometrika. 2006;93:491-507.
#'
#' @export
adjust_TSBH <- function(rawp, alpha = 0.05) {
  m <- length(rawp)
  index <- order(rawp)
  spval <- rawp[index]
  tmp <- spval

  # Step 1: BH adjusted p-values
  for (i in (m - 1):1)
    tmp[i] <- min(tmp[i + 1], min((m / i) * spval[i], 1, na.rm = TRUE),
                  na.rm = TRUE)

  # Step 2: estimate m0 from first-pass rejections at alpha/(1+alpha)
  h0.TSBH <- m - sum(tmp < alpha / (1 + alpha), na.rm = TRUE)
  tmp <- pmin(tmp * h0.TSBH / m, 1)

  result <- numeric(m)
  result[index] <- tmp
  return(result)
}
