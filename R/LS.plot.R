#' Plot Lomb–Scargle periodogram for a single feature
#'
#' @param object Seurat object
#' @param time.col Numeric metadata column name containing sampling times
#' @param feature Single gene/feature name to plot
#' @param assay Assay name (default "RNA")
#' @param slot Slot name (default "data")
#' @param center Logical; if TRUE, subtract the mean after windowing
#' @param window.func String name of window function (e.g. "hanning") or NULL
#' @param f.min Numeric; minimum frequency
#' @param f.max Numeric; maximum frequency
#' @param n.bins Integer; number of frequency bins
#' @return A ggplot object of frequency vs. power
#' @import signal
#' @importFrom ggplot2 ggplot aes geom_line labs theme_bw
#' @export
LS.plot <- function(object, time.col, feature,
                    assay = "RNA", slot = "data", center = FALSE,
                    window.func = NULL, f.min = 0, f.max = 2.0,
                    n.bins = 500) {
  if (length(feature) != 1) stop("`feature` must be a single gene name.")

  np <- reticulate::import("numpy")
  ats <- reticulate::import("astropy.timeseries")

  md <- object@meta.data
  if (!(time.col %in% colnames(md))) stop(paste("Metadata column", time.col, "not found."))
  time.vec <- md[[time.col]]
  expr.mat <- Seurat::GetAssayData(object, assay = assay, slot = slot)
  if (!(feature %in% rownames(expr.mat))) {
    stop(sprintf("Feature '%s' not found in assay '%s'.", feature, assay))
  }
  expr.vec <- as.numeric(expr.mat[feature, ])
  keep <- is.finite(time.vec) & is.finite(expr.vec)
  tvec <- time.vec[keep]; svec <- expr.vec[keep]
  ord <- order(tvec); tvec <- tvec[ord]; svec <- svec[ord]
  N <- length(tvec)
  if (!is.null(window.func)) { wf <- match.fun(window.func); wv <- wf(N) } else { wv <- rep(1, N) }
  ywin <- svec * wv; if (center) ywin <- ywin - mean(ywin)
  t.py <- np$array(tvec); y.py <- np$array(ywin)
  freqs <- np$linspace(f.min, f.max, as.integer(n.bins))
  ls.obj <- ats$LombScargle(t.py, y.py)
  power <- ls.obj$power(freqs)
  df <- data.frame(Frequency = py_to_r(freqs), Power = py_to_r(power))
  ggplot(df, aes(x = Frequency, y = Power)) +
    geom_line() + labs(title = feature, x = "Frequency", y = "Power") +
    theme_bw()
}
