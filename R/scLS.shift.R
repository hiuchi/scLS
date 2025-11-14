#' Perform dual Lomb–Scargle and compute distance between power spectra
#' for separate group1 and group2 Seurat objects
#'
#' This function assumes that group1 and group2 cells are stored in two separate Seurat
#' objects (group1.object and group2.object). Each object has its own pseudotime column
#' in metadata. For each feature (gene), Lomb–Scargle periodograms are computed
#' separately for group1 and group2 along their respective pseudotime axes, and a
#' distance between the two power spectra is calculated. A null distribution of
#' distances is obtained by permuting pseudotimes within group1 and group2 separately,
#' and a Gaussian-based p-value is reported.
#'
#' @param group1.object A Seurat object for the group1 group, with pseudotime in metadata
#'   and expression in an assay.
#' @param group2.object A Seurat object for the group2 group, with pseudotime in metadata
#'   and expression in an assay.
#' @param time.col1 Character. Metadata column name in \code{group1.object} containing
#'   pseudotime values for group1.
#' @param time.col2 Character. Metadata column name in \code{group2.object} containing
#'   pseudotime values for group2.
#' @param features Optional. Vector of feature (gene) names to analyze. If NULL
#'   (default), the intersection of \code{Seurat::VariableFeatures(group1.object)} and
#'   \code{Seurat::VariableFeatures(group2.object)} is used, further intersected with
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
#'       (group1 vs group2).}
#'     \item{p}{Gaussian-based p-value from the fitted null distribution
#'       (NA if computation fails).}
#'   }
#'
#' @import signal
#' @importFrom stats dist pnorm sd
#' @importFrom tibble tibble
#' @export
scLS.shift <- function(
    group1.object,
    group2.object,
    time.col1,
    time.col2,
    features    = NULL,
    assay       = "RNA",
    slot        = "data",
    center      = FALSE,
    window.func = NULL,
    f.min       = 0,
    f.max       = 2.0,
    n.bins      = 500,
    dist.method = "canberra",
    n.perm      = 100,
    n.cores     = 1,
    seed        = 8
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
  meta_group1 <- group1.object@meta.data
  meta_group2 <- group2.object@meta.data

  if (!time.col1 %in% colnames(meta_group1)) {
    stop("time.col1 must exist in group1.object@meta.data.")
  }
  if (!time.col2 %in% colnames(meta_group2)) {
    stop("time.col2 must exist in group2.object@meta.data.")
  }

  group1_cells <- colnames(group1.object)
  group2_cells <- colnames(group2.object)

  t1_vec <- as.numeric(meta_group1[group1_cells, time.col1])
  t2_vec <- as.numeric(meta_group2[group2_cells, time.col2])

  group1_valid <- !is.na(t1_vec)
  group2_valid <- !is.na(t2_vec)

  group1_cells <- group1_cells[group1_valid]
  group2_cells <- group2_cells[group2_valid]
  t1_vec       <- t1_vec[group1_valid]
  t2_vec       <- t2_vec[group2_valid]

  ## ---- Expression matrices & features ----
  group1_expr <- Seurat::GetAssayData(group1.object, assay = assay, slot = slot)
  group2_expr <- Seurat::GetAssayData(group2.object, assay = assay, slot = slot)

  common_genes <- intersect(rownames(group1_expr), rownames(group2_expr))

  if (is.null(features)) {
    vf_group1 <- Seurat::VariableFeatures(group1.object)
    vf_group2 <- Seurat::VariableFeatures(group2.object)
    features  <- intersect(vf_group1, vf_group2)
  }
  features <- intersect(features, common_genes)

  if (length(features) < 1) {
    stop("No common features found between group1.object and group2.object (after intersecting VariableFeatures and rownames).")
  }

  ## Restrict expression matrices to valid cells
  group1_expr <- group1_expr[features, group1_cells, drop = FALSE]
  group2_expr <- group2_expr[features, group2_cells, drop = FALSE]

  ## Common frequency grid
  freqs_py <- np$linspace(f.min, f.max, as.integer(n.bins))

  ## ---- Per-feature computation ----
  compute_feature <- function(gene) {
    ## Expression vectors
    sig1 <- as.numeric(group1_expr[gene, ])
    sig2 <- as.numeric(group2_expr[gene, ])

    ## Variance-based handling of (quasi-)constant groups
    eps_var  <- 1e-8   # threshold to regard a signal as (almost) constant
    eps_jitt <- 1e-3   # noise scale added to constant signals

    v1 <- stats::var(sig1, na.rm = TRUE)
    v2 <- stats::var(sig2, na.rm = TRUE)

    ## Both groups almost constant -> no dynamic information; return NA
    if (!is.na(v1) && v1 < eps_var && !is.na(v2) && v2 < eps_var) {
      return(tibble::tibble(
        gene     = gene,
        distance = NA_real_,
        p        = NA_real_
      ))
    }

    ## Only one group almost constant -> add small noise to stabilise LS
    if (!is.na(v1) && v1 < eps_var) {
      sig1 <- sig1 + stats::rnorm(length(sig1), mean = 0, sd = eps_jitt)
    }
    if (!is.na(v2) && v2 < eps_var) {
      sig2 <- sig2 + stats::rnorm(length(sig2), mean = 0, sd = eps_jitt)
    }

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

    ## Lomb–Scargle for group1 and group2
    ls1 <- ats$LombScargle(t1_py, y1_py)
    ls2 <- ats$LombScargle(t2_py, y2_py)
    p1  <- reticulate::py_to_r(ls1$power(freqs_py))
    p2  <- reticulate::py_to_r(ls2$power(freqs_py))

    ## Guard against non-finite powers
    if (!all(is.finite(p1)) || !all(is.finite(p2))) {
      return(tibble::tibble(
        gene     = gene,
        distance = NA_real_,
        p        = NA_real_
      ))
    }

    obs_dist <- as.numeric(stats::dist(rbind(p1, p2), method = dist.method))
    if (!is.finite(obs_dist)) {
      obs_dist <- NA_real_
    }

    ## Permutation null: shuffle pseudotimes within group1 and group2 separately
    perm_d <- replicate(n.perm, {
      t1p_py <- np$array(sample(t1_vec))
      t2p_py <- np$array(sample(t2_vec))

      d1 <- reticulate::py_to_r(ats$LombScargle(t1p_py, y1_py)$power(freqs_py))
      d2 <- reticulate::py_to_r(ats$LombScargle(t2p_py, y2_py)$power(freqs_py))

      if (!all(is.finite(d1)) || !all(is.finite(d2))) {
        return(NA_real_)
      }
      as.numeric(stats::dist(rbind(d1, d2), method = dist.method))
    })

    mu0 <- mean(perm_d, na.rm = TRUE)
    sd0 <- stats::sd(perm_d, na.rm = TRUE)

    if (is.na(sd0) || sd0 == 0 || is.na(mu0)) {
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
