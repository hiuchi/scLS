#' Perform dual Lomb–Scargle and compute distance between power spectra
#' for separate WT and KO Seurat objects
#'
#' This function assumes that WT and KO cells are stored in two separate Seurat
#' objects (wt.object and ko.object). Each object has its own pseudotime column
#' in metadata. For each feature (gene), Lomb–Scargle periodograms are computed
#' separately for WT and KO along their respective pseudotime axes, and a
#' distance between the two power spectra is calculated. A null distribution of
#' distances is obtained by permuting pseudotimes within WT and KO separately,
#' and a Gaussian-based p-value is reported.
#'
#' @param wt.object A Seurat object for the WT group, with pseudotime in metadata
#'   and expression in an assay.
#' @param ko.object A Seurat object for the KO group, with pseudotime in metadata
#'   and expression in an assay.
#' @param time.col1 Character. Metadata column name in \code{wt.object} containing
#'   pseudotime values for WT.
#' @param time.col2 Character. Metadata column name in \code{ko.object} containing
#'   pseudotime values for KO.
#' @param features Optional. Vector of feature (gene) names to analyze. If NULL
#'   (default), the intersection of \code{Seurat::VariableFeatures(wt.object)} and
#'   \code{Seurat::VariableFeatures(ko.object)} is used, further intersected with
#'   genes present in both expression matrices.
#' @param assay Character. Assay name to extract expression data from for both
#'   objects. Default: "RNA".
#' @param slot Character. Slot in assay to pull data from for both objects.
#'   Default: "data".
#' @param center Logical; if TRUE subtract the mean after windowing (for each
#'   group separately). Default: FALSE.
#' @param window.func Character, function, or NULL; name of window function in
#'   'signal' package (e.g., "hanning"), or a custom function taking \code{n}
#'   and returning a length-\code{n} vector, or NULL for no windowing.
#'   Default: NULL.
#' @param f.min Numeric; minimum frequency. Default: 0.
#' @param f.max Numeric; maximum frequency. Default: 2.0.
#' @param n.bins Integer; number of frequency bins. Default: 500.
#' @param dist.method Character; distance method for \code{stats::dist()}.
#'   Default: "canberra".
#' @param n.perm Integer; number of permutations for null distribution
#'   (pseudotime permutation within each group). Default: 100.
#' @param seed Integer or NULL; random seed for reproducibility of permutations
#'   (R and NumPy). Default: 8.
#' @param n.cores Integer; number of cores for parallel execution across
#'   features (uses \code{parallel::mclapply}). Default: 1.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{gene}{Feature name analyzed.}
#'     \item{distance}{Observed distance between the two power spectra
#'       (WT vs KO).}
#'     \item{p}{Gaussian-based p-value from the fitted null distribution
#'       (NA if computation fails).}
#'   }
#'
#' @import signal
#' @importFrom stats dist pnorm sd
#' @importFrom tibble tibble
#' @export
scLS.shift <- function(
  wt.object,
  ko.object,
  time.col1,
  time.col2,
  features   = NULL,
  assay      = "RNA",
  slot       = "data",
  center     = FALSE,
  window.func = NULL,
  f.min      = 0,
  f.max      = 2.0,
  n.bins     = 500,
  dist.method = "canberra",
  n.perm     = 100,
  n.cores    = 1,
  seed       = 8
) {
  ## ---- Dependency checks ----
  if (!requireNamespace("signal", quietly = TRUE)) {
    stop("Package 'signal' is required.")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.")
  }
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required.")
  }
  if (n.cores > 1 && !requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required for n.cores > 1.")
  }

  ## ---- RNG & Python setup ----
  if (!is.null(seed)) {
    set.seed(seed)
  }
  np  <- reticulate::import("numpy")
  if (!is.null(seed) && !is.null(np$random)) {
    np$random$seed(as.integer(seed))
  }
  np$seterr(divide = "ignore", invalid = "ignore")
  ats <- reticulate::import("astropy.timeseries")

  ## ---- Metadata & pseudotime ----
  meta_wt <- wt.object@meta.data
  meta_ko <- ko.object@meta.data

  if (!time.col1 %in% colnames(meta_wt)) {
    stop("time.col1 must exist in wt.object@meta.data.")
  }
  if (!time.col2 %in% colnames(meta_ko)) {
    stop("time.col2 must exist in ko.object@meta.data.")
  }

  wt_cells <- colnames(wt.object)
  ko_cells <- colnames(ko.object)

  t1_vec <- as.numeric(meta_wt[wt_cells, time.col1])
  t2_vec <- as.numeric(meta_ko[ko_cells, time.col2])

  wt_valid <- !is.na(t1_vec)
  ko_valid <- !is.na(t2_vec)

  wt_cells <- wt_cells[wt_valid]
  ko_cells <- ko_cells[ko_valid]
  t1_vec   <- t1_vec[wt_valid]
  t2_vec   <- t2_vec[ko_valid]

  ## ---- Expression matrices & features ----
  wt_expr <- Seurat::GetAssayData(wt.object, assay = assay, slot = slot)
  ko_expr <- Seurat::GetAssayData(ko.object, assay = assay, slot = slot)

  common_genes <- intersect(rownames(wt_expr), rownames(ko_expr))

  if (is.null(features)) {
    vf_wt <- Seurat::VariableFeatures(wt.object)
    vf_ko <- Seurat::VariableFeatures(ko.object)
    features <- intersect(vf_wt, vf_ko)
  }
  features <- intersect(features, common_genes)

  if (length(features) < 1) {
    stop("No common features found between wt.object and ko.object (after intersecting VariableFeatures and rownames).")
  }

  ## Restrict expression matrices to valid cells
  wt_expr <- wt_expr[features, wt_cells, drop = FALSE]
  ko_expr <- ko_expr[features, ko_cells, drop = FALSE]

  ## Common frequency grid
  freqs_py <- np$linspace(f.min, f.max, as.integer(n.bins))

  ## ---- Per-feature computation ----
  compute_feature <- function(gene) {
    ## Expression vectors
    sig1 <- as.numeric(wt_expr[gene, ])
    sig2 <- as.numeric(ko_expr[gene, ])

    n1 <- length(sig1)
    n2 <- length(sig2)

    ## Windowing
    if (is.character(window.func)) {
      if (!exists(window.func, where = asNamespace("signal"))) {
        stop(paste0("Window function '", window.func, "' not found in 'signal' package."))
      }
      w_fun <- get(window.func, envir = asNamespace("signal"))
      w1 <- w_fun(n1)
      w2 <- w_fun(n2)
    } else if (is.function(window.func)) {
      w1 <- window.func(n1)
      w2 <- window.func(n2)
    } else {
      w1 <- rep(1, n1)
      w2 <- rep(1, n2)
    }

    y1 <- sig1 * w1
    y2 <- sig2 * w2
    if (center) {
      y1 <- y1 - mean(y1)
      y2 <- y2 - mean(y2)
    }

    ## Convert to numpy arrays
    t1_py <- np$array(t1_vec)
    t2_py <- np$array(t2_vec)
    y1_py <- np$array(y1)
    y2_py <- np$array(y2)

    ## Lomb–Scargle for WT and KO
    ls1 <- ats$LombScargle(t1_py, y1_py)
    ls2 <- ats$LombScargle(t2_py, y2_py)
    p1  <- reticulate::py_to_r(ls1$power(freqs_py))
    p2  <- reticulate::py_to_r(ls2$power(freqs_py))

    obs_dist <- as.numeric(stats::dist(rbind(p1, p2), method = dist.method))

    ## Permutation null: shuffle pseudotimes within WT and KO separately
    perm_d <- replicate(n.perm, {
      t1p_py <- np$array(sample(t1_vec))
      t2p_py <- np$array(sample(t2_vec))

      d1 <- reticulate::py_to_r(ats$LombScargle(t1p_py, y1_py)$power(freqs_py))
      d2 <- reticulate::py_to_r(ats$LombScargle(t2p_py, y2_py)$power(freqs_py))

      as.numeric(stats::dist(rbind(d1, d2), method = dist.method))
    })

    mu0 <- mean(perm_d)
    sd0 <- stats::sd(perm_d)

    if (is.na(sd0) || sd0 == 0) {
      if (is.na(obs_dist)) {
        pval <- NA_real_
      } else {
        pval <- if (obs_dist > mu0) 0 else 1
      }
    } else {
      ## Right-sided p-value: P(D >= obs_dist)
      pval <- 1 - stats::pnorm(obs_dist, mean = mu0, sd = sd0)
    }

    tibble::tibble(gene = gene, distance = obs_dist, p = pval)
  }

  ## ---- Parallel wrapper ----
  if (n.cores > 1) {
    chunks <- split(features, cut(seq_along(features), breaks = n.cores, labels = FALSE))
    res_list <- parallel::mclapply(
      chunks,
      function(chunk) lapply(chunk, compute_feature),
      mc.cores = n.cores
    )
    results <- do.call(c, res_list)
  } else {
    results <- lapply(features, compute_feature)
  }

  do.call(rbind, results)
}
