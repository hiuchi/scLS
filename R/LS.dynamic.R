#' Perform Lomb–Scargle analysis with Astropy
#'
#' @param object Seurat object
#' @param time.col Numeric metadata column name containing sampling times
#' @param feature Gene/feature name to analyze; if NULL, use VariableFeatures(object)
#' @param assay Assay name (default "RNA")
#' @param slot Slot name (default "data")
#' @param center Logical; if TRUE, subtract the mean after windowing
#' @param window.func String name of window function (e.g. "hanning") or NULL
#' @param f.min Numeric; minimum frequency
#' @param f.max Numeric; maximum frequency
#' @param n.bins Integer; number of frequency bins
#' @param fap.method Character; method for Astropy FAP estimation (default "baluev")
#' @param n.cores Integer; number of cores for parallel execution (default 1)
#' @return A tibble with columns Feature, PeakFrequency, and PeakFAP
#' @import signal
#' @importFrom Seurat GetAssayData VariableFeatures
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
#' @importFrom stats setNames
#' @importFrom signal hanning
#' @importFrom reticulate import py_to_r
#' @export
LS.dynamic <- function(object, time.col, feature = NULL, assay = "RNA", slot = "data",
                       center = FALSE, window.func = NULL, f.min = 0, f.max = 2.0,
                       n.bins = 500, fap.method = "baluev", n.cores = 1) {
  if (!requireNamespace("signal", quietly = TRUE)) stop("Package 'signal' is required.")
  if (!requireNamespace("parallel", quietly = TRUE)) stop("Package 'parallel' is required.")
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object.")
  if (is.null(feature)) feature <- VariableFeatures(object)

  np <- reticulate::import("numpy")
  ats <- reticulate::import("astropy.timeseries")

  md <- object@meta.data
  if (!(time.col %in% colnames(md))) stop(paste("Metadata column", time.col, "not found."))
  time.vec <- md[[time.col]]
  expr.mat <- Seurat::GetAssayData(object, assay = assay, slot = slot)
  process_feature <- function(feat) {
    removed <- setNames(0L, feat)
    if (!(feat %in% rownames(expr.mat))) {
      warning(paste("Feature", feat, "not found; skipped."))
      return(list(result = NULL, removed = removed))
    }
    sig.vec <- as.numeric(expr.mat[feat, ])
    keep <- is.finite(time.vec) & is.finite(sig.vec)
    removed[] <- sum(!keep)
    if (!any(keep)) {
      warning(paste("Feature", feat, "has no finite observations after filtering; skipped."))
      return(list(result = NULL, removed = removed))
    }
    tvec <- time.vec[keep]; svec <- sig.vec[keep]
    ord <- order(tvec); tvec <- tvec[ord]; svec <- svec[ord]
    N <- length(tvec)
    if (!is.null(window.func)) {
      wf <- match.fun(window.func); wv <- wf(N)
    } else {
      wv <- rep(1, N)
    }
    ywin <- svec * wv
    if (center) ywin <- ywin - mean(ywin)
    t.py <- np$array(tvec); y.py <- np$array(ywin)
    freqs <- np$linspace(f.min, f.max, as.integer(n.bins))
    ls.obj <- ats$LombScargle(t.py, y.py)
    power <- ls.obj$power(freqs)
    fap.vec <- ls.obj$false_alarm_probability(power, method = fap.method)
    freq.r <- py_to_r(freqs); fap.r <- py_to_r(fap.vec)
    idx <- which.max(power)
    tib <- tibble::tibble(
      Feature = feat,
      PeakFrequency = freq.r[idx],
      PeakFAP = fap.r[idx]
    )
    list(result = tib, removed = removed)
  }

  if (n.cores > 1) {
    res_list <- parallel::mclapply(feature, process_feature, mc.cores = n.cores)
  } else {
    res_list <- lapply(feature, process_feature)
  }

  removed.list <- lapply(res_list, function(x) x$removed)
  if (length(removed.list) > 0) {
    removed.counts <- unlist(removed.list, use.names = TRUE)
  } else {
    removed.counts <- setNames(integer(0), character(0))
  }
  res_flat <- lapply(res_list, function(x) x$result)
  res_flat <- Filter(Negate(is.null), res_flat)

  total_removed <- sum(removed.counts)
  n_feats <- sum(removed.counts > 0)
  if (total_removed > 0) message(paste(total_removed, "NA/Inf removed across", n_feats, "features"))

  if (!length(res_flat)) {
    return(tibble::tibble(Feature = character(), PeakFrequency = numeric(), PeakFAP = numeric()))
  }

  results <- dplyr::bind_rows(res_flat)
  return(results)
}
