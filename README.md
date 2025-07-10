## Table of Contents

1. **[Introduction & Environment Setup](#introduction--environment-setup)**  
   *Motivation, hardware / OS assumptions, Conda-based Python environment, and R-package installation*

2. **[Dataset Acquisition & Pseudotime Computation](#dataset-acquisition--pseudotime-computation)**  
   *Download PBMC example dataset, load into Seurat, run Slingshot, add pseudotime metadata, and perform quick QC plots*

3. **[Lomb–Scargle Analysis & Spectral Visualisation](#lombscargle-analysis--spectral-visualisation)**  
   *Single-gene and genome-wide `LS.dynamic` runs, parallelisation options, and `LS.plot` power-spectrum figures*

4. **[Comparing Two Trajectories with `LS.shift`](#comparing-two-trajectories-with-lsshift)**  
   *Distance metrics, permutation *p*-values, multi-gene analyses, and CPU-scaling considerations*
