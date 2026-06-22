#' Compare multiple testing correction methods
#'
#' Applies BH, BY, ABH, TSBH, and CORA corrections to the same set of
#' p-values and returns a summary table.
#'
#' @param rawp Numeric vector of raw p-values.
#' @param cor_matrix Square matrix of absolute pairwise correlations.
#'   Required for CORA. If NULL, CORA is skipped.
#' @param alpha Significance level (default 0.05).
#'
#' @return A data.frame with columns: method, n_rejected, pct_of_BH.
#'
#' @examples
#' set.seed(42)
#' m <- 500
#' rawp <- c(runif(50, 0, 0.005), runif(450, 0, 1))
#' cor_mat <- diag(m)
#' # Add block correlation among first 50 genes
#' cor_mat[1:50, 1:50] <- 0.3
#' diag(cor_mat) <- 1
#'
#' cora_compare(rawp, cor_mat)
#'
#' @export
cora_compare <- function(rawp, cor_matrix = NULL, alpha = 0.05) {

  methods <- list(
    BH   = adjust_BH(rawp),
    BY   = adjust_BY(rawp),
    ABH  = adjust_ABH(rawp),
    TSBH = adjust_TSBH(rawp, alpha)
  )

  if (!is.null(cor_matrix)) {
    methods$CORA <- adjust_CORA(rawp, cor_matrix)
  }

  n_BH <- sum(methods$BH <= alpha)

  result <- data.frame(
    method     = names(methods),
    n_rejected = sapply(methods, function(adj) sum(adj <= alpha)),
    pct_of_BH  = sapply(methods, function(adj) {
      if (n_BH > 0) round(sum(adj <= alpha) / n_BH * 100, 1) else NA
    }),
    stringsAsFactors = FALSE
  )
  rownames(result) <- NULL

  return(result)
}


#' Summary of CORA adjustment
#'
#' Provides a detailed summary of the CORA adjustment including the number
#' of rejections, reduction vs BH, and top genes.
#'
#' @param rawp Numeric vector of raw p-values (named for gene identification).
#' @param cor_matrix Square matrix of absolute pairwise correlations.
#' @param alpha Significance level (default 0.05).
#' @param top_n Number of top genes to display (default 20).
#'
#' @return A list with components:
#'   \item{comparison}{Data frame comparing all methods}
#'   \item{reduction}{Percentage reduction in DEGs: CORA vs BH}
#'   \item{top_genes}{Data frame of top genes with raw and adjusted p-values}
#'   \item{overlap}{Overlap between CORA and BH top-N gene sets}
#'
#' @examples
#' set.seed(42)
#' m <- 200
#' rawp <- c(runif(30, 0, 0.005), runif(170, 0, 1))
#' names(rawp) <- paste0("gene", 1:m)
#' cor_mat <- diag(m)
#' cor_mat[1:30, 1:30] <- 0.4
#' diag(cor_mat) <- 1
#'
#' res <- cora_summary(rawp, cor_mat)
#' print(res$comparison)
#' cat("DEG reduction:", res$reduction, "%\n")
#'
#' @export
cora_summary <- function(rawp, cor_matrix, alpha = 0.05, top_n = 20) {

  adj_bh   <- adjust_BH(rawp)
  adj_cora <- adjust_CORA(rawp, cor_matrix)

  n_bh   <- sum(adj_bh   <= alpha)
  n_cora <- sum(adj_cora <= alpha)

  # Reduction
  reduction <- if (n_bh > 0) round((1 - n_cora / n_bh) * 100, 1) else 0

  # Comparison table
  comparison <- cora_compare(rawp, cor_matrix, alpha)

  # Top genes
  ord <- order(rawp)
  top_idx <- ord[1:min(top_n, length(rawp))]

  top_genes <- data.frame(
    gene      = if (!is.null(names(rawp))) names(rawp)[top_idx] else top_idx,
    raw_p     = rawp[top_idx],
    adj_BH    = adj_bh[top_idx],
    adj_CORA  = adj_cora[top_idx],
    sig_BH    = adj_bh[top_idx]   <= alpha,
    sig_CORA  = adj_cora[top_idx] <= alpha,
    stringsAsFactors = FALSE
  )

  # Overlap at different top-N
  overlaps <- data.frame()
  for (n in c(50, 100, 200)) {
    if (n > length(rawp)) next
    top_bh   <- ord[1:n]
    top_cora <- order(adj_cora)[1:n]
    ov <- length(intersect(top_bh, top_cora))
    overlaps <- rbind(overlaps, data.frame(
      top_n   = n,
      overlap = ov,
      pct     = round(ov / n * 100, 1)
    ))
  }

  return(list(
    comparison = comparison,
    reduction  = reduction,
    top_genes  = top_genes,
    overlap    = overlaps
  ))
}
