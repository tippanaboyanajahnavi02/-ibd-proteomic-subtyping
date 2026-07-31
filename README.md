# ProSubIBD: Proteomic Profiling to Identify and Validate IBD Subtypes

MSc Health Data Science Thesis
University of Birmingham, 2026
Student: Tippanaboyana Jahnavi (2989054)
Supervisor: Dr Animesh Acharjee

## Overview

Ensemble deep learning consensus
clustering pipeline for proteomic
subtyping of inflammatory bowel
disease (IBD) across two independent
international cohorts.

## Key Results

- Two robust proteomic subtypes
  identified across all four cohorts
- Hyperinflammatory subtype:
  TNFα/NFκB, IL-6/JAK/STAT3,
  IFN-γ signalling elevated
- Silhouette scores: 0.777--0.846
- Random Forest AUC: 0.983--0.998
- External validation in IBDome
  (Germany): r = 0.735, Δ = 0.014

## Repository Structure

| Notebook | Description |
|---|---|
| 01_preprocessing | QC, MICE imputation, INT normalisation |
| 02_traditional_clustering | Baseline clustering methods |
| 03_autoencoders | AE, DAE, VAE, β-VAE, Batch-VAE training |
| 04_dec | Deep Embedded Clustering fine-tuning |
| 05_consensus | Weighted co-occurrence consensus |
| 06_biological_characterisation | DE, RF, stability analysis |
| 07_limma | R-based limma differential abundance |
| 08_ORA | MSigDB Hallmark pathway enrichment |
| 09_PPI_network_analysis | STRING PPI network construction |
| ibdome_cd_pipeline | IBDome Crohn's Disease pipeline |
| ibdome_uc_pipeline | IBDome Ulcerative Colitis pipeline |
| Figures | Main thesis figures |
| Supplementary_figures | Supplementary figures S1-S12 |

## Software Requirements

- Python 3.12.3
- PyTorch 2.12.0
- scikit-learn 1.4.2
- R 4.3.1 with limma 3.58.1

Full software versions in
Supplementary Table S11.

## Data Access

- UK Biobank: application required
  at ukbiobank.ac.uk
- IBDome: publicly available at
  ibdome.de

## Citation

Tippanaboyana Jahnavi (2026).
ProSubIBD: Proteomic Profiling to
Identify and Validate IBD Subtypes.
MSc Health Data Science thesis,
University of Birmingham.
