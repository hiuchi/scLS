# Interpreting scLS p-values in human hematopoiesis

This vignette demonstrates how to use `scLS.dynamic()` to rank
pseudotime-associated genes in a well-characterized human hematopoietic
differentiation dataset. We use the human CD34+ bone marrow scRNA-seq data
analyzed by Setty et al. (2019), which includes Palantir pseudotime and branch
probabilities.

Genes are ranked by `PeakFAP`, the p-value returned by `scLS.dynamic()`. Genes
with smaller values are prioritized for trend visualization and comparison with
known lineage markers.

## Dataset

Setty et al. applied Palantir to human CD34+ bone marrow cells and provide
processed `scanpy` AnnData objects for three biological replicates:

- [Replicate 1](https://s3.amazonaws.com/dp-lab-data-public/palantir/human_cd34_bm_rep1.h5ad)
- [Replicate 2](https://s3.amazonaws.com/dp-lab-data-public/palantir/human_cd34_bm_rep2.h5ad)
- [Replicate 3](https://s3.amazonaws.com/dp-lab-data-public/palantir/human_cd34_bm_rep3.h5ad)

The objects include log-normalized expression values, Palantir pseudotime,
branch probabilities, and low-dimensional coordinates. The example below uses
Replicate 1 for a compact, reproducible demonstration.

Relevant references:

- Palantir GitHub repository: <https://github.com/dpeerlab/Palantir>
- Palantir manuscript: <https://doi.org/10.1038/s41587-019-0068-4>

## 1. Install dependencies

Install the R packages used in this example:

```r
install.packages(c("anndata", "curl", "dplyr", "ggplot2", "patchwork", "remotes"))
remotes::install_github("hiuchi/scLS")
```

`scLS` uses Python's Astropy through `reticulate`. Create a Python environment
that contains Astropy before running `scLS.dynamic()`:

```bash
conda create -n scLS-env python=3.12 astropy anndata -c conda-forge
conda activate scLS-env
```

## 2. Load the Palantir AnnData object

This code downloads Replicate 1 to an absolute temporary path and converts it to
a Seurat object. The Seurat object is created from the raw count matrix stored in
the AnnData raw layer, and normalization is performed with Seurat.

Run the following code in a fresh R session. `reticulate` must be attached to
`scLS-env` before any package initializes Python.

```r
library(reticulate)
use_condaenv("scLS-env", required = TRUE)

library(anndata)
library(curl)
library(Matrix)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scLS)

rep1_url <- "https://s3.amazonaws.com/dp-lab-data-public/palantir/human_cd34_bm_rep1.h5ad"
rep1_file <- file.path(tempdir(), "human_cd34_bm_rep1.h5ad")

curl_download(rep1_url, rep1_file)
adata <- read_h5ad(rep1_file)

raw_counts <- as(Matrix::t(adata$raw$X), "CsparseMatrix")
rownames(raw_counts) <- rownames(adata$raw$var)
colnames(raw_counts) <- rownames(adata$obs)

human_cd34 <- CreateSeuratObject(
  counts = raw_counts,
  meta.data = adata$obs,
  project = "Setty_CD34_BM_Rep1"
)
```

Add Palantir t-SNE coordinates and branch probabilities to the Seurat object:

```r
tsne <- adata$obsm[["tsne"]]
rownames(tsne) <- rownames(adata$obs)
colnames(tsne) <- c("tSNE_1", "tSNE_2")

human_cd34[["tsne"]] <- CreateDimReducObject(
  embeddings = tsne,
  key = "tSNE_",
  assay = "RNA"
)

branch_probs <- as.data.frame(adata$obsm[["palantir_branch_probs"]])
colnames(branch_probs) <- make.names(adata$uns[["palantir_branch_probs_cell_types"]])
rownames(branch_probs) <- rownames(adata$obs)

human_cd34 <- AddMetaData(human_cd34, metadata = branch_probs)
```

Normalize the raw counts and select variable genes for the transcriptome-wide
`scLS.dynamic()` scan using the standard Seurat workflow.

```r
human_cd34 <- NormalizeData(
  human_cd34,
  verbose = FALSE
)

human_cd34 <- FindVariableFeatures(
  human_cd34,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)
```

## 3. Select one lineage for a clear dynamic-expression example

`scLS.dynamic()` tests for expression structure along one supplied pseudotime
axis. For a branching differentiation system, apply it to a biologically
interpretable lineage or trajectory segment when the goal is lineage-specific
interpretation.

Here we use the branch-probability column corresponding to the erythroid branch.
The exact column name is taken from the Palantir metadata.

```r
erythroid_col <- colnames(branch_probs)[grepl("ery", colnames(branch_probs), ignore.case = TRUE)][1]
erythroid_cells <- colnames(human_cd34)[human_cd34[[erythroid_col]][, 1] > 0.5]

erythroid_obj <- subset(human_cd34, cells = erythroid_cells)

pt <- erythroid_obj$palantir_pseudotime
erythroid_obj$pseudotime_scLS <- (pt - min(pt, na.rm = TRUE)) /
  (max(pt, na.rm = TRUE) - min(pt, na.rm = TRUE))
```

This branch-restricted analysis keeps the interpretation simple: genes with
small `PeakFAP` values are genes whose expression varies along the erythroid
pseudotime ordering.

## 4. Run scLS.dynamic()

Start with a small marker panel to confirm the interpretation. Known erythroid
genes such as `GATA1`, `KLF1`, `TAL1`, `HBA1`, `HBA2`, `HBB`, `GYPA`, and
`ALAS2` are useful genes to inspect for this lineage.

```r
erythroid_markers <- c("GATA1", "KLF1", "TAL1", "GYPA", "HBA1", "HBA2", "HBB", "ALAS2")
erythroid_markers <- intersect(erythroid_markers, rownames(erythroid_obj))

marker_result <- scLS.dynamic(
  object = erythroid_obj,
  time.col = "pseudotime_scLS",
  feature = erythroid_markers,
  assay = "RNA",
  slot = "data",
  center = TRUE,
  window.func = "hanning",
  f.min = 0.01,
  f.max = 2,
  n.bins = 500,
  fap.method = "baluev"
)

marker_result %>%
  mutate(q_value = p.adjust(PeakFAP, method = "BH")) %>%
  arrange(PeakFAP)
```

Then scan the variable genes and rank genes by `PeakFAP`:

```r
all_result <- scLS.dynamic(
  object = erythroid_obj,
  time.col = "pseudotime_scLS",
  assay = "RNA",
  slot = "data",
  center = TRUE,
  window.func = "hanning",
  f.min = 0.01,
  f.max = 2,
  n.bins = 500,
  fap.method = "baluev",
  n.cores = 4
)

ranked_result <- all_result %>%
  mutate(q_value = p.adjust(PeakFAP, method = "BH")) %>%
  arrange(PeakFAP)

ranked_result %>%
  select(Feature, PeakFrequency, PeakFAP, q_value) %>%
  head(20)
```

In a local run using Replicate 1 with an erythroid branch-probability threshold
of 0.5, erythroid-associated genes such as `BLVRB`, `KLF1`, `GATA1`, `AHSP`,
`ANK1`, `TFR2`, `HBB`, and `HBA1` were highly ranked. This is the intended use
of the p-value ranking: it prioritizes genes for marker comparison and
trajectory-level visualization.

## 5. Plot representative genes along pseudotime

The next step is visual inspection. A low `PeakFAP` value should lead to a
clear, interpretable expression pattern along pseudotime.

```r
plot_gene_trend <- function(object, gene) {
  expr <- as.numeric(GetAssayData(object, assay = "RNA", slot = "data")[gene, ])

  data.frame(
    pseudotime = object$pseudotime_scLS,
    expression = expr
  ) %>%
    ggplot(aes(x = pseudotime, y = expression)) +
    geom_point(size = 0.4, alpha = 0.35) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.8) +
    labs(title = gene, x = "Palantir pseudotime", y = "Expression") +
    theme_bw()
}

genes_to_plot <- intersect(
  c("GATA1", "KLF1", "AHSP", "HBB", "BLVRB", "ANK1", "TFR2"),
  rownames(erythroid_obj)
)
trend_plots <- lapply(genes_to_plot, plot_gene_trend, object = erythroid_obj)

wrap_plots(trend_plots)
```

For this dataset, top-ranked genes should be interpreted by checking whether
they recover known hematopoietic differentiation markers or transcriptional
programs. In an erythroid branch, useful genes to inspect include erythroid
regulators and markers such as `GATA1`, `KLF1`, `AHSP`, `ANK1`, `HBB`, and
`HBA1`. For a myeloid branch, a similar workflow can be applied after selecting
the corresponding branch-probability column, where genes such as `MPO`, `ELANE`,
`LYZ`, `SPI1`, `S100A8`, and `S100A9` may be useful markers to inspect.

## 6. Practical interpretation

Use the `scLS.dynamic()` result as a ranking table:

- Sort genes by `PeakFAP` or by a multiple-testing adjusted value such as
  `p.adjust(PeakFAP, method = "BH")`.
- Inspect the top-ranked gene trends along pseudotime.
- Compare top genes with known markers, transcription factors, and pathway-level
  annotations.
- Treat the result as a first-pass screening step. The biological conclusion
  should come from the combination of ranking, visualization, lineage context,
  and external biological knowledge.

For branching datasets such as hematopoiesis, avoid interpreting a single global
run across all branches as branch-specific evidence. When the question is
lineage-specific, first select a lineage or trajectory segment, then run
`scLS.dynamic()` on that subset.

## 7. Why this example is useful

Human hematopoietic differentiation is a well-studied system. That makes it a
good practical check for the end-user utility of `scLS.dynamic()` p-values:
instead of asking whether the method produces a different numerical ranking from
another tool, this workflow asks whether the ranking helps users recover and
inspect biologically meaningful trajectory-associated genes in a familiar
differentiation system.
