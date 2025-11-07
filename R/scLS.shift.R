#' Perform dual Lomb–Scargle and compute distance between power spectra for Seurat objects
#'
#' This function takes a Seurat object containing two pseudotime metadata columns and expression data
#' for features (genes). It computes Lomb–Scargle periodograms along each pseudotime axis
#' for each feature and returns a tibble with feature names, distances between power spectra, and Gaussian-based p-values from permutations.
#'
#' @param object A Seurat object with pseudotime in metadata and expression in an assay.
#' @param time.col1 Character. Name of metadata column containing the first pseudotime values.
#' @param time.col2 Character. Name of metadata column containing the second pseudotime values.
#' @param features Optional. Vector of feature (gene) names to analyze. If NULL (default), all Seurat::VariableFeatures(object) are used.
#' @param assay Character. Assay name to extract expression data from. Default: "RNA".
#' @param slot Character. Slot in assay to pull data from. Default: "data".
#' @param center Logical; if TRUE subtract the mean after windowing. Default: FALSE.
#' @param window.func Character or function or NULL; name of window function in 'signal' package (e.g., "hanning"), or a custom function, or NULL. Default: NULL.
#' @param f.min Numeric; minimum frequency. Default: 0.
#' @param f.max Numeric; maximum frequency. Default: 2.0.
#' @param n.bins Integer; number of frequency bins. Default: 500.
#' @param dist.method Character; distance method for dist(). Default: "canberra".
#' @param n.perm Integer; number of permutations for null distribution. Default: 100.
#' @param seed Integer or NULL; random seed for reproducibility of permutations. Default: 8.
#' @param n.cores Integer; number of cores for parallel execution. Default: 1.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{gene}{Feature name analyzed.}
#'     \item{distance}{Observed distance between the two power spectra.}
#'     \item{p}{Gaussian-based p-value from fitted null distribution (NA if computation fails).}
#'   }
#'
#' @import signal
#' @importFrom stats dist pnorm sd
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot aes geom_line labs theme_bw
#' @export
scLS.shift <- function(object, time.col1, time.col2, features = NULL, assay = "RNA", slot = "data",
                     center = FALSE, window.func = NULL, f.min = 0, f.max = 2.0,
                     n.bins = 500, dist.method = "canberra", n.perm = 100, n.cores = 1, seed = 8) {
  if (!requireNamespace("signal", quietly = TRUE)) stop("Package 'signal' is required.")
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required.")
  if (!requireNamespace("reticulate", quietly = TRUE)) stop("Package 'reticulate' is required.")
  if (!requireNamespace("parallel", quietly = TRUE)) stop("Package 'parallel' is required for n.cores > 1.")

  set.seed(seed)
  np <- reticulate::import("numpy")
  if (!is.null(seed) && !is.null(np$random)) np$random$seed(as.integer(seed))
  np$seterr(divide = "ignore", invalid = "ignore")
  ats <- reticulate::import("astropy.timeseries")

  meta <- object@meta.data
  if (!(time.col1 %in% colnames(meta)) || !(time.col2 %in% colnames(meta))) {
    stop("time.col1 and time.col2 must exist in metadata.")
  }
  cells_all <- colnames(object)
  time1_all <- meta[cells_all, time.col1]
  time2_all <- meta[cells_all, time.col2]
  valid_cells <- cells_all[!is.na(time1_all) & !is.na(time2_all)]
  if (length(valid_cells) < 3) stop("Not enough cells with valid pseudotimes.")

  if (is.null(features)) {
    features <- Seurat::VariableFeatures(object)
    if (length(features) < 1) stop("No variable features found.")
  }

  expr_mat <- Seurat::GetAssayData(object, assay = assay, slot = slot)

  compute_feature <- function(gene) {
    sig_all <- as.numeric(expr_mat[gene, cells_all])
    sig <- sig_all[match(valid_cells, cells_all)]
    n <- length(sig)

    if (is.character(window.func)) {
      if (!exists(window.func, where = asNamespace("signal"))) stop(paste0("Window function '", window.func, "' not found in 'signal' package."))
      w <- get(window.func, envir = asNamespace("signal"))(n)
    } else if (is.function(window.func)) {
      w <- window.func(n)
    } else {
      w <- rep(1, n)
    }
    y <- sig * w
    if (center) y <- y - mean(y)

    t1_py <- np$array(meta[valid_cells, time.col1])
    t2_py <- np$array(meta[valid_cells, time.col2])
    y_py <- np$array(y)
    freqs_py <- np$linspace(f.min, f.max, as.integer(n.bins))

    ls1 <- ats$LombScargle(t1_py, y_py)
    ls2 <- ats$LombScargle(t2_py, y_py)
    p1 <- reticulate::py_to_r(ls1$power(freqs_py))
    p2 <- reticulate::py_to_r(ls2$power(freqs_py))
    obs_dist <- as.numeric(dist(rbind(p1, p2), method = dist.method))

    perm_d <- replicate(n.perm, {
      t1p <- np$array(sample(meta[valid_cells, time.col1]))
      t2p <- np$array(sample(meta[valid_cells, time.col2]))
      d1 <- reticulate::py_to_r(ats$LombScargle(t1p, y_py)$power(freqs_py))
      d2 <- reticulate::py_to_r(ats$LombScargle(t2p, y_py)$power(freqs_py))
      as.numeric(dist(rbind(d1, d2), method = dist.method))
    })
    mu0 <- mean(perm_d)
    sd0 <- sd(perm_d)

    if (is.na(sd0) || sd0 == 0) {
      pval <- if (is.na(obs_dist)) NA_real_ else if (obs_dist > mu0) 0 else 1
    } else {
      pval <- 1 - pnorm(obs_dist, mu0, sd0)
    }

    tibble::tibble(gene = gene, distance = obs_dist, p = pval)
  }

  if (n.cores > 1) {
    chunks <- split(features, cut(seq_along(features), breaks = n.cores, labels = FALSE))
    res_list <- parallel::mclapply(chunks, function(chunk) lapply(chunk, compute_feature), mc.cores = n.cores)
    results <- do.call(c, res_list)
  } else {
    results <- lapply(features, compute_feature)
  }

  do.call(rbind, results)
}
