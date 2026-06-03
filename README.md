# CORA: COrrelation-Redundancy-Aware FDR Adjustment

[![R package](https://img.shields.io/badge/R-package-blue.svg)](https://github.com/zjanna/CORA)

## Overview

**CORA** (COrrelation-Redundancy-Aware) modifies the Benjamini-Hochberg FDR
control procedure by incorporating pairwise correlations into the
rejection threshold. Hypotheses strongly correlated with already-rejected ones
receive a stricter threshold, reducing redundancy in the list of
discoveries. While designed and benchmarked on genomic data (microarray
and RNA-seq), CORA is applicable to any multiple testing problem where a
correlation matrix among test statistics is available.

## Installation

```r
# Install CORA from GitHub
devtools::install_github("zjanna/CORA")
```

## Quick Start

```r
library(CORA)

# Your p-values and correlation matrix
rawp <- compute_your_pvalues(data)       # e.g., from t-tests or DESeq2
cor_matrix <- abs(cor(t(expression_data))) # absolute gene-gene correlations

# CORA-adjusted p-values
adj_cora <- adjust_CORA(rawp, cor_matrix)

# Compare all methods
cora_compare(rawp, cor_matrix, alpha = 0.05)
#>   method n_rejected pct_of_BH
#> 1     BH       1756     100.0
#> 2     BY        403      22.9
#> 3    ABH       2182     124.3
#> 4   TSBH       2030     115.6
#> 5   CORA       1519      86.5

# Detailed summary
res <- cora_summary(rawp, cor_matrix)
cat("DEG reduction vs BH:", res$reduction, "%\n")
```

## How it works

For each hypothesis at rank *i* in the p-value ranking, CORA computes:

$$c_i = \sum_{k=1}^{i-1} |r_{\sigma(k), \sigma(i)}|$$

the sum of absolute correlations with all previously ranked hypotheses.
The adjusted p-value becomes:

$$\tilde{p}_{(i)}^{\text{CORA}} = \min_{j \geq i} \left\{\frac{m}{j - c_j} \cdot p_{(j)}\right\}$$

**Key properties:**
- When correlations are zero: CORA = BH (identical results)
- When correlations are positive: CORA is more conservative than BH
- CORA rejections are always a subset of BH rejections (proven)
- Empirically validated: FDR_CORA <= alpha in all tested scenarios (10 simulation scenarios)

## When to use CORA

CORA is most valuable when:
- Test statistics are organized in correlated modules (e.g., co-expression pathways)
- The number of discoveries is large (> 200)
- You want to reduce redundancy in the discovery list
- You are selecting biomarker panels (top-N genes)

## Citation

If you use CORA in your research, please cite:

> Zyprych-Walczak J. CORA: a COrrelation-Redundancy-Aware FDR adjustment
> with genomic applications. BMC Bioinformatics. 2026. (submitted)

## License

GPL (>= 3)
