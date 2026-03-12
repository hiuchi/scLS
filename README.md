# scLS

`scLS` applies Lomb-Scargle analysis to single-cell trajectories stored in Seurat
objects. This README is written as a short vignette: it shows how to install the
package, connect R to a Python environment that provides Astropy, compute
pseudotime with Slingshot, run `scLS.dynamic()` on one trajectory, and run
`scLS.shift()` to compare two trajectories.

The workflow below uses an example PBMC dataset so you can run the full analysis
end to end before switching to your own data.

## What the two main functions do

- `scLS.dynamic()`: analyze one pseudotime axis and find genes with periodic or
  oscillatory structure along that trajectory.
- `scLS.shift()`: compare two trajectories and quantify how different each gene's
  Lomb-Scargle spectrum is between them.

## 1. Install scLS and its dependencies

`scLS` uses Python's Astropy through `reticulate`, so you need both R packages
and a Python environment.

### 1.1 Create the Python environment and download the example dataset

Run these commands in a terminal.

```bash
# Create an isolated Python 3.12 environment
conda create -n scLS-env python=3.12 -c conda-forge

# Activate the environment
conda activate scLS-env

# Install Python dependencies used through reticulate
conda install astropy -c conda-forge

# Download the example Seurat object into the current directory
curl -L -o pbmc10k_mye_small_velocyto.rds \
  https://zenodo.org/records/10944066/files/pbmc10k_mye_small_velocyto.rds
```

### 1.2 Install the R packages

From this point onward, run commands in R.

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("hiuchi/scLS")
devtools::install_github("huayc09/SeuratExtend")
```

### 1.3 Bind R to the Python environment

```r
library(reticulate)

use_condaenv("scLS-env", required = TRUE)

# Optional sanity check
py_config()
reticulate::import("astropy.timeseries")
```

## 2. Load the example object and compute pseudotime

If your Seurat object already contains numeric pseudotime columns in
`object@meta.data`, you can skip to Section 3.

```r
library(Seurat)
library(SeuratExtend)
library(tidyverse)
library(scLS)

# Read the example object
mye_small <- readRDS("pbmc10k_mye_small_velocyto.rds")

# Compute Slingshot pseudotime
mye_small <- RunSlingshot(
  Seu = mye_small,
  group.by = "cluster",
  start.clus = "Mono CD14"
)

# Copy the Slingshot pseudotime matrix into meta.data
sling <- mye_small@misc$slingshot$PCA$SlingPseudotime
mye_small[[]] <- cbind(mye_small[[]], sling)

# The new metadata columns are now available
colnames(mye_small[[]])[grepl("^slingPseudotime_", colnames(mye_small[[]]))]
```

Optional QC:

```r
DimPlot2(mye_small, features = colnames(sling), cols = "C")
GeneTrendCurve.Slingshot(mye_small, features = c("CD14", "FCGR3A"))
```

## 3. Run `scLS.dynamic()` on one trajectory

`scLS.dynamic()` takes one Seurat object and one pseudotime column. Start with a
single gene to confirm that the pipeline is working, then scale up to a larger
feature set.

### 3.1 Single-gene example

```r
result_cd14 <- scLS.dynamic(
  object = mye_small,
  time.col = "slingPseudotime_1",
  feature = "CD14",
  center = TRUE,
  window.func = "hanning",
  f.min = 0,
  f.max = 2,
  n.bins = 500,
  fap.method = "baluev"
)

result_cd14
```

The returned tibble contains:

- `Feature`: gene name
- `PeakFrequency`: frequency with maximum Lomb-Scargle power
- `PeakFAP`: false alarm probability at that peak

### 3.2 Scan many genes

If `feature = NULL`, `scLS.dynamic()` uses `VariableFeatures(object)`. That is a
good default for a first transcriptome-wide scan.

```r
result_all <- scLS.dynamic(
  object = mye_small,
  time.col = "slingPseudotime_1",
  center = TRUE,
  window.func = "hanning",
  f.min = 0,
  f.max = 2,
  n.bins = 500,
  fap.method = "baluev",
  n.cores = 4
)

result_all %>%
  arrange(PeakFAP) %>%
  head(10)
```

Use this table to identify genes that show strong structure along one
trajectory. In practice, you will usually sort by `PeakFAP` and then inspect the
top genes in more detail.

## 4. Run `scLS.shift()` to compare two trajectories

`scLS.shift()` expects two Seurat objects, one for each trajectory being
compared. A simple pattern is to start from one Seurat object that contains
multiple Slingshot lineages and then create lineage-specific objects by keeping
cells with finite pseudotime in each lineage.

### 4.1 Create lineage-specific Seurat objects

```r
lineage1_obj <- subset(
  mye_small,
  cells = colnames(mye_small)[is.finite(mye_small$slingPseudotime_1)]
)

lineage2_obj <- subset(
  mye_small,
  cells = colnames(mye_small)[is.finite(mye_small$slingPseudotime_2)]
)
```

### 4.2 Compare one gene between the two trajectories

```r
shift_cd14 <- scLS.shift(
  group1.object = lineage1_obj,
  group2.object = lineage2_obj,
  time.col1 = "slingPseudotime_1",
  time.col2 = "slingPseudotime_2",
  features = "CD14",
  center = TRUE,
  window.func = "hanning",
  n.perm = 20
)

shift_cd14
```

The returned tibble contains:

- `gene`: gene name
- `distance`: distance between the two Lomb-Scargle spectra
- `p`: permutation-based significance estimate for that distance

### 4.3 Compare a small gene panel first

`scLS.shift()` is more computationally expensive than `scLS.dynamic()` because it
builds a permutation-based null distribution. For a README-sized example, start
with a small panel of variable genes.

```r
shift_features <- intersect(
  VariableFeatures(lineage1_obj),
  VariableFeatures(lineage2_obj)
)[1:50]

shift_panel <- scLS.shift(
  group1.object = lineage1_obj,
  group2.object = lineage2_obj,
  time.col1 = "slingPseudotime_1",
  time.col2 = "slingPseudotime_2",
  features = shift_features,
  center = TRUE,
  window.func = "hanning",
  n.perm = 20,
  n.cores = 4
)

shift_panel %>%
  arrange(p) %>%
  head(10)
```

Genes with larger `distance` and smaller `p` are stronger candidates for
trajectory-specific dynamics. For a larger production run, remove the
`features = shift_features` line and increase `n.perm`.

## 5. How to adapt this to your own data

To use `scLS` on another dataset, you need:

1. A Seurat object with expression values in the assay you want to analyze.
2. One or more numeric pseudotime columns stored in `object@meta.data`.
3. For `scLS.shift()`, two Seurat objects representing the two trajectories you
   want to compare.

Typical substitutions are:

- replace `mye_small` with your own Seurat object
- replace `slingPseudotime_1` and `slingPseudotime_2` with your own pseudotime
  column names
- replace `"CD14"` with genes of interest
- adjust `f.min`, `f.max`, `n.bins`, and `n.cores` to suit your dataset

If you already computed pseudotime with another method, you can still use
`scLS.dynamic()` and `scLS.shift()` as long as the pseudotime values are numeric
and stored in metadata.

## 6. Minimal working summary

If you want the shortest possible path:

1. Install `scLS`, `SeuratExtend`, and Astropy.
2. Load a Seurat object and add a numeric pseudotime column.
3. Run `scLS.dynamic()` on one pseudotime axis.
4. Split cells into two trajectory-specific Seurat objects.
5. Run `scLS.shift()` to compare the two trajectories.

That is the core `scLS` workflow.
