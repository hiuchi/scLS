### 1. Environment Setup <a name="introduction--environment-setup"></a>

#### 1.1 Create, activate, and populate the Conda environment *(Terminal)*

```bash
# 1) Create an isolated Python 3.12 environment
conda create -n scLS-env python=3.12 -c conda-forge

# 2) Initialise Conda for zsh and reload the shell
conda init zsh
exec $SHELL        # …or simply restart your Terminal

# 3) Activate the new environment
conda activate scLS-env

# 4) Install core Python dependency (Astropy)
conda install astropy -c conda-forge

# 5) Download the example PBMC dataset into the current directory
curl -L -o pbmc10k_mye_small_velocyto.rds \
  https://zenodo.org/records/10944066/files/pbmc10k_mye_small_velocyto.rds

```

#### 1.2 Install the required R packages
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("hiuchi/scLS")          # scLS (Lomb–Scargle utilities)
devtools::install_github("huayc09/SeuratExtend")  # SeuratExtend (Slingshot helpers)
```

#### 1.3 Bind R to the Conda environment
```
library(reticulate)

use_condaenv("scLS-env", required = TRUE)

# Verify that Python packages are available
np  <- reticulate::import("numpy")
ats <- reticulate::import("astropy.timeseries")
```

### 2. Dataset Acquisition & Pseudotime Computation <a name="dataset-acquisition--pseudotime-computation"></a>

> **All commands below run in _R_.**

```r
# Load libraries
library(Seurat)
library(SeuratExtend)
library(tidyverse)

# 2.1  Read the example PBMC object downloaded in Step 1.1
mye_small <- readRDS("pbmc10k_mye_small_velocyto.rds")

# 2.2  Compute pseudotime with Slingshot
mye_small <- RunSlingshot(
  object    = mye_small,
  group.by  = "cluster",        # clustering column in meta.data
  start.clus = "Mono CD14"      # root cluster
)

# 2.3  Move Slingshot pseudotime matrix into meta.data
sling <- mye_small@misc$slingshot$PCA$SlingPseudotime
mye_small[[]] <- cbind(mye_small[[]], sling)

# 2.4  Quick QC: visualise pseudotime embeddings
DimPlot2(mye_small, features = colnames(sling), cols = "C")

# 2.5  Optional exploratory gene-trend plots
GeneTrendCurve.Slingshot(mye_small, features = c("CD14", "FCGR3A"))

GeneTrendHeatmap.Slingshot(
  mye_small,
  features = c("CD14", VariableFeatures(mye_small)[1:10]),
  lineage  = "slingPseudotime_2"
)
```

### 3. Lomb–Scargle Analysis & Spectral Visualisation <a name="lombscargle-analysis--spectral-visualisation"></a>

```r
### 3.1  Single-gene Lomb–Scargle periodogram (CD14)
result_cd14 <- scLS.dynamic(
  object      = mye_small,
  time.col    = "slingPseudotime_1",
  feature     = "CD14",
  center      = TRUE,
  window.func = hanning,
  f.min       = 0,
  f.max       = 2,
  n.bins      = 500,
  fap.method  = "baluev"
)

### 3.2  Genome-wide scan on variable features (parallel, 4 cores)
result_all <- scLS.dynamic(
  object      = mye_small,
  time.col    = "slingPseudotime_1",
  center      = TRUE,
  window.func = hanning,
  f.min       = 0,
  f.max       = 2,
  n.bins      = 500,
  fap.method  = "baluev",
  n.cores     = 4
)

### 3.3  Visualise the power spectrum for CD14
g <- scLS.plot(
  object      = mye_small,
  time.col    = "slingPseudotime_1",
  feature     = "CD14",
  center      = TRUE,
  window.func = "hanning",
  f.min       = 0,
  f.max       = 1,
  n.bins      = 100
)
print(g)   # displays ggplot object
```

### 4. Comparing Two Trajectories with LS.shift <a name="comparing-two-trajectories-with-lsshift"></a>
```r
### 4.1  Distance / permutation p-value for a single gene
shift_cd14 <- scLS.shift(
  object      = mye_small,
  time.col1   = "slingPseudotime_1",
  time.col2   = "slingPseudotime_2",
  feature     = "CD14",
  center      = TRUE,
  window.func = "hanning"
)

# 4.2  Transcriptome-wide comparison (12 cores)
shift_all <- scLS.shift(
  object      = mye_small,
  time.col1   = "slingPseudotime_1",
  time.col2   = "slingPseudotime_2",
  center      = TRUE,
  window.func = "hanning",
  n.cores     = 12
)
```
