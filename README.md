# ProSubIBD: Proteomic Profiling to Identify and Validate IBD Subtypes

**MSc Health Data Science Thesis**
University of Birmingham, 2026
Student: Tippanaboyana Jahnavi (ID: 2989054)
Supervisor: Dr Animesh Acharjee
Institute: Department of Cancer and Genomic Sciences

---

## Project Overview

ProSubIBD is an ensemble deep learning
consensus clustering framework for
unsupervised proteomic subtyping of
inflammatory bowel disease (IBD).
The pipeline integrates five autoencoder
architectures with Deep Embedded
Clustering and three clustering
algorithms to generate 15 independent
solutions, combined into a
silhouette-weighted consensus.
Two robust proteomic subtypes were
identified and validated across two
independent international cohorts:
a Hyperinflammatory subtype
characterised by TNF-α/NF-κB,
IL-6/JAK/STAT3, and IFN-γ signalling,
and a Quiescent subtype with markedly
lower inflammatory burden.

### Key Results

| Cohort | Disease | n | Proteins | Silhouette | Stability | AUC |
|--------|---------|---|----------|-----------|-----------|-----|
| UK Biobank | CD | 215 | 991 | 0.846 | 0.893 | 0.983 |
| UK Biobank | UC | 430 | 991 | 0.777 | 0.844 | 0.985 |
| IBDome | CD | 201 | 61 | 0.834 | 0.886 | 0.992 |
| IBDome | UC | 132 | 61 | 0.827 | 0.882 | 0.998 |

- Zero unstably assigned patients
  in three of four cohorts
- 100% agreement between Random
  Forest and limma in all cohorts
- Cross-national validation:
  r = 0.735 vs reference r = 0.749
  (Δ = 0.014)

---
Apply at: https://www.ukbiobank.ac.uk
Platform: Olink Explore 3072
2,923 proteins measured via
Proximity Extension Assay
Patients: Identified using ICD-10
K50 = Crohn's disease
K51 = Ulcerative colitis
Access: Data available upon approved
application and data access
agreement
Note: UK Biobank data cannot be
shared directly. Researchers
must apply independently.

### IBDome
Website: https://www.ibdome.de
Access: Publicly available
Centres: Charité Berlin +
University Hospital Erlangen
Germany
Citation: Plattner et al. (2025)
bioRxiv doi:10.1101/
2025.03.26.645544
Download: ibdome_cd_preprocessed.csv
ibdome_uc_preprocessed.csv
61 proteins, 333 patients

---

## Computational Environment

### Hardware — BlueBEAR HPC
Cluster: BlueBEAR HPC
Institution: University of Birmingham
Edgbaston, Birmingham
United Kingdom

GPU Node Specification:
GPU: NVIDIA A100-SXM4 40GB
GPU Memory: 40 GB HBM2e
GPU Count: 1 per job
CUDA: 13.0

CPU Specification:
Processor: AMD EPYC 7742
Cores: 8 per job (requested)
RAM: 64 GB per job (requested)

Storage:
Home: /rds/homes/[username]/
Scratch: /rds/projects/[project]/
Quota: 1 TB

Scheduler: SLURM Workload Manager
Queue: bbgpu (GPU queue)
Max runtime: 24 hours per job


### Software Environment

Operating System: Linux (Ubuntu 22.04)
Python: 3.12.3
GCCcore-13.3.0
CUDA: 13.0
cuDNN: 8.9.7
PyTorch: 2.12.0+cu130
R: 4.3.1


---

## Complete Setup Instructions

### Step 1 — Access BlueBEAR

```bash
# SSH into BlueBEAR from your terminal
ssh [your_username]@bluebear.bham.ac.uk

# Enter your university password
# when prompted

# Navigate to your home directory
cd /rds/homes/[your_username]/
```

### Step 2 — Load Required Modules

```bash
# Always run this before anything else
module purge

# Load core modules
module load bluebear
module load Python/3.12.3-GCCcore-13.3.0
module load CUDA/13.0.0
module load GCC/13.3.0

# Verify modules loaded correctly
module list
```

### Step 3 — Clone Repository

```bash
# Clone this repository
git clone https://github.com/\
tippanaboyanajahnavi02/\
-ibd-proteomic-subtyping.git

# Enter the project directory
cd -ibd-proteomic-subtyping

# Create directories for outputs
mkdir -p data tables figures logs
```

### Step 4 — Create Python Environment

```bash
# Create virtual environment
python -m venv prosubIBD_env

# Activate it
source prosubIBD_env/bin/activate

# Verify Python version
python --version
# Expected: Python 3.12.3

# Upgrade pip
pip install --upgrade pip setuptools wheel
```

### Step 5 — Install All Python Packages

```bash
# Install all required packages
# with exact versions for reproducibility

pip install \
    numpy==2.4.6 \
    pandas==3.0.3 \
    scipy==1.17.1 \
    scikit-learn==1.8.0 \
    torch==2.12.0 \
    torchvision \
    umap-learn==0.5.12 \
    hdbscan==0.8.38 \
    networkx==3.6.1 \
    gseapy==1.2.1 \
    matplotlib==3.10.9 \
    seaborn==0.13.2 \
    statsmodels==0.14.6 \
    ipykernel==6.29.5 \
    jupyter \
    jupyterlab \
    tqdm \
    joblib \
    requests \
    py2cytoscape \
    plotly

# Verify key packages installed
python -c "
import torch
import numpy as np
import sklearn
print('PyTorch:', torch.__version__)
print('CUDA available:',
      torch.cuda.is_available())
print('NumPy:', np.__version__)
print('sklearn:', sklearn.__version__)
"
```

### Step 6 — Install R and Packages

```bash
# Load R module
module load R/4.3.1-foss-2023a

# Open R
R

# Inside R console install packages:
install.packages(
    "BiocManager",
    repos="https://cran.r-project.org")

BiocManager::install("limma")

install.packages(
    c("mice",
      "dplyr",
      "tidyr",
      "ggplot2",
      "readr"),
    repos="https://cran.r-project.org")

# Verify installations
library(limma)
library(mice)
packageVersion("limma")
# Expected: 3.58.1
packageVersion("mice")
# Expected: 3.16.0

# Exit R
q()
```

### Step 7 — Register Jupyter Kernel

```bash
# Register virtual environment
# as a Jupyter kernel
python -m ipykernel install \
    --user \
    --name prosubIBD \
    --display-name "ProSubIBD (Python 3.12)"

# Verify kernel registered
jupyter kernelspec list
```

### Step 8 — Request Interactive GPU Session

```bash
# Request an interactive GPU node
# for running notebooks directly
srun \
    --account=[your_account_code] \
    --qos=bbgpu \
    --gres=gpu:a100:1 \
    --time=08:00:00 \
    --mem=64G \
    --cpus-per-task=8 \
    --nodes=1 \
    --pty bash

# Once on GPU node load modules again
module purge
module load bluebear
module load Python/3.12.3-GCCcore-13.3.0
module load CUDA/13.0.0

# Activate environment
source /rds/homes/[username]/\
prosubIBD_env/bin/activate

# Verify GPU is visible
python -c "
import torch
print('GPU available:',
      torch.cuda.is_available())
print('GPU name:',
      torch.cuda.get_device_name(0))
print('GPU memory:',
      torch.cuda.get_device_properties(0)\
      .total_memory / 1e9, 'GB')
"
# Expected output:
# GPU available: True
# GPU name: NVIDIA A100-SXM4-40GB
# GPU memory: 40.54 GB
```

### Step 9 — Launch Jupyter Notebook

```bash
# On the GPU node start Jupyter
jupyter notebook \
    --no-browser \
    --port=8888 \
    --ip=$(hostname -i)

# Note the URL shown e.g.:
# http://10.xx.xx.xx:8888/?token=xxxx

# On your LOCAL machine open a
# NEW terminal and run SSH tunnel:
ssh -L 8888:[gpu_node_hostname]:8888 \
    [username]@bluebear.bham.ac.uk

# Then open in your browser:
# http://localhost:8888
# Enter the token from above
```

### Step 10 — SLURM Batch Job Script

For long training runs without
interactive session use this
batch script. Save as
`run_prosubIBD.sh`:

```bash
#!/bin/bash
#SBATCH --account=[your_account_code]
#SBATCH --qos=bbgpu
#SBATCH --gres=gpu:a100:1
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --job-name=prosubIBD_train
#SBATCH --output=logs/train_%j.out
#SBATCH --error=logs/train_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=[your_email]

# Load modules
module purge
module load bluebear
module load Python/3.12.3-GCCcore-13.3.0
module load CUDA/13.0.0

# Activate environment
source /rds/homes/${USER}/\
prosubIBD_env/bin/activate

# Change to project directory
cd /rds/homes/${USER}/\
-ibd-proteomic-subtyping

echo "Job started: $(date)"
echo "Running on: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name \
    --format=csv,noheader)"

# Execute notebooks in order
for nb in \
    01_preprocessing.ipynb \
    02_traditional_clustering.ipynb \
    03_autoencoders.ipynb \
    04_dec.ipynb \
    05_consensus.ipynb \
    06_biological_characterisation.ipynb \
    08_ORA.ipynb \
    09_PPI_network_analysis.ipynb
do
    echo "Running: $nb"
    jupyter nbconvert \
        --to notebook \
        --execute \
        --inplace \
        --ExecutePreprocessor.timeout=7200 \
        "$nb"
    echo "Completed: $nb"
done

echo "Job finished: $(date)"
```

Submit the job:

```bash
sbatch run_prosubIBD.sh

# Monitor job status
squeue -u [your_username]

# View output log
tail -f logs/train_[jobid].out

# Cancel job if needed
scancel [jobid]
```

---

## Running the Analysis

### Notebook Execution Order

Run ALL notebooks in this exact order:

STEP 1 — Preprocessing
Notebook: 01_preprocessing.ipynb
Runtime: ~15 minutes (CPU)
Input: Raw UK Biobank data
Output: X_cd_int.npy (215x991)
X_uc_int.npy (430x991)
meta_cd_final.csv
meta_uc_final.csv

STEP 2 — Baseline Clustering
Notebook: 02_traditional_clustering.ipynb
Runtime: ~10 minutes (CPU)
Input: X_cd_int.npy, X_uc_int.npy
Output: Baseline comparison metrics

STEP 3 — Autoencoder Training
Notebook: 03_autoencoders.ipynb
Runtime: ~90 minutes (A100 GPU)
~8 hours (CPU only)
Input: X_cd_int.npy, X_uc_int.npy
Output: Latent representations
for all 5 architectures

STEP 4 — DEC Fine-tuning
Notebook: 04_dec.ipynb
Runtime: ~30 minutes (A100 GPU)
~3 hours (CPU only)
Input: Latent representations
Output: DEC-refined latents
latent_cd_dae_dec.npy
latent_uc_dae_dec.npy

STEP 5 — Consensus Clustering
Notebook: 05_consensus.ipynb
Runtime: ~20 minutes (CPU)
Input: All latent representations
Output: consensus_matrix_cd.npy
consensus_matrix_uc.npy
labels_cd_consensus.npy
labels_uc_consensus.npy
stability_cd.npy
stability_uc.npy

STEP 6 — Biological Characterisation
Notebook: 06_biological_characterisation.ipynb
Runtime: ~30 minutes (CPU)
Input: Consensus labels + matrices
Output: RF importance tables
Silhouette profiles

STEP 7 — Limma DE Analysis (R)
Notebook: 07_limma.ipynb
Runtime: ~10 minutes (R kernel)
Input: X_cd_int.npy, labels
Output: cd_limma_M3.csv (991x7)
uc_limma_M3.csv (991x7)
NOTE: Select R kernel in Jupyter
not Python kernel

STEP 8 — Pathway Enrichment
Notebook: 08_ORA.ipynb
Runtime: ~15 minutes (CPU)
Requires internet connection
for Enrichr API
Input: DE protein lists
Output: ORA Hallmark tables

STEP 9 — PPI Network Analysis
Notebook: 09_PPI_network_analysis.ipynb
Runtime: ~20 minutes (CPU)
Requires internet for STRING
Input: DE protein lists
Output: PPI network figures

STEP 10 — IBDome CD Validation
Notebook: ibdome_cd_pipeline.ipynb
Runtime: ~45 minutes (A100 GPU)
Input: IBDome CD data
Output: ibdome_cd_validation_v2.json
r = 0.735, Δ = 0.014

STEP 11 — IBDome UC Validation
Notebook: ibdome_uc_pipeline.ipynb
Runtime: ~30 minutes (A100 GPU)
Input: IBDome UC data
Output: ibdome_uc_validation_v2.json
r = 0.689, Δ = 0.060

STEP 12 — Generate Main Figures
Notebook: Figures.ipynb
Runtime: ~15 minutes (CPU)
Output: figure1 through figure5

STEP 13 — Generate Supplementary Figures
Notebook: Supplementary_figures.ipynb
Runtime: ~20 minutes (CPU)
Output: figS1 through figS12


---

## Reproducibility

### Random Seeds

Every stochastic operation uses
**seed = 42**. Add this at the top
of every notebook:

```python
import random
import numpy as np
import torch
import os

SEED = 42

# Python
random.seed(SEED)

# NumPy
np.random.seed(SEED)

# PyTorch CPU
torch.manual_seed(SEED)

# PyTorch GPU
torch.cuda.manual_seed(SEED)
torch.cuda.manual_seed_all(SEED)

# CuDNN determinism
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

# OS level
os.environ['PYTHONHASHSEED'] = str(SEED)

print(f'All seeds set to {SEED}')
```

### Device Configuration

```python
# Standard device setup used in
# all GPU training notebooks

import torch

DEVICE = torch.device(
    'cuda' if torch.cuda.is_available()
    else 'cpu')

print(f'Using device: {DEVICE}')

if DEVICE.type == 'cuda':
    print(f'GPU: '
          f'{torch.cuda.get_device_name(0)}')
    print(f'Memory: '
          f'{torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB')
    # Expected on BlueBEAR:
    # GPU: NVIDIA A100-SXM4-40GB
    # Memory: 40.5 GB
```

### Expected Runtimes

On NVIDIA A100-SXM4 40GB (BlueBEAR):
Autoencoder pretraining: ~90 min
DEC fine-tuning: ~30 min
Consensus clustering: ~20 min
Biological analysis: ~45 min
IBDome validation: ~75 min
Figure generation: ~35 min
TOTAL: ~5 hours

On CPU only (no GPU):
Autoencoder pretraining: ~8-12 hours
All other steps: same as above
TOTAL: ~15 hours


### Important Reproducibility Notes
GPU non-determinism
Even with seed=42, minor numerical
differences may occur across
different GPU hardware due to
floating point parallelism.
Results should agree to 3+
decimal places.
Pre-computed outputs
All .npy output files are saved.
You can skip GPU training and
load directly: Z = np.load('data/latent_cd_dae_dec.npy')
labels = np.load('data/labels_cd_consensus.npy')
Internet required
Notebooks 08 (ORA via Enrichr API)
and 09 (STRING PPI database)
require internet connection.
BlueBEAR nodes have internet access.
R kernel for notebook 07
Notebook 07_limma.ipynb must be
run using the R kernel not Python.
In Jupyter: Kernel > Change Kernel

R 4.3.1

Data paths
Update DATA_DIR in each notebook
to match your file locations:
DATA_DIR = '/rds/homes/[username]/data/'

---

## Software Versions (Full List)

| Software | Version | Purpose |
|----------|---------|---------|
| Python | 3.12.3 | Main language |
| PyTorch | 2.12.0+cu130 | Deep learning |
| CUDA | 13.0 | GPU acceleration |
| cuDNN | 8.9.7 | GPU neural nets |
| numpy | 2.4.6 | Array operations |
| pandas | 3.0.3 | Data handling |
| scipy | 1.17.1 | Statistics |
| scikit-learn | 1.8.0 | ML algorithms |
| umap-learn | 0.5.12 | UMAP projection |
| hdbscan | 0.8.38 | Clustering |
| networkx | 3.6.1 | PPI networks |
| gseapy | 1.2.1 | Pathway analysis |
| matplotlib | 3.10.9 | Visualisation |
| seaborn | 0.13.2 | Statistical plots |
| statsmodels | 0.14.6 | Statistics |
| ipykernel | 6.29.5 | Jupyter kernel |
| R | 4.3.1 | Statistical computing |
| limma (R) | 3.58.1 | Differential expression |
| mice (R) | 3.16.0 | MICE imputation |
| BiocManager (R) | 1.30.23 | Bioconductor |
| STRING database | v12.0 | PPI interactions |
| MSigDB | v2023.2 | Gene sets |
| GCC | 13.3.0 | C compiler |
| SLURM | 22.05 | Job scheduler |

---

## Output Files

Running the complete pipeline
generates the following files.
All files are saved automatically
by each notebook.

---

### data/

Core preprocessed matrices,
cluster labels, latent
representations, and metadata
for all four cohorts.

#### UK Biobank — Preprocessed Matrices

| File | Shape | Description |
|------|-------|-------------|
| `X_cd_int.npy` | 215 × 991 | CD protein matrix after INT normalisation |
| `X_uc_int.npy` | 430 × 991 | UC protein matrix after INT normalisation |
| `meta_cd_final.csv` | 215 × 12 | CD patient metadata with cluster assignments |
| `meta_uc_final.csv` | 430 × 12 | UC patient metadata with cluster assignments |

#### UK Biobank — Cluster Labels

| File | Shape | Description |
|------|-------|-------------|
| `labels_cd_consensus.npy` | (215,) | Final consensus cluster labels CD |
| `labels_uc_consensus.npy` | (430,) | Final consensus cluster labels UC |
| `labels_cd_gmm_best.npy` | (215,) | Independent GMM labels for ARI validation CD (ARI = 0.963) |
| `labels_uc_gmm_best.npy` | (430,) | Independent GMM labels for ARI validation UC (ARI = 0.972) |

#### UK Biobank — Consensus Matrices

| File | Shape | Description |
|------|-------|-------------|
| `consensus_matrix_cd.npy` | 215 × 215 | Silhouette-weighted co-occurrence matrix CD |
| `consensus_matrix_uc.npy` | 430 × 430 | Silhouette-weighted co-occurrence matrix UC |

#### UK Biobank — Latent Representations

| File | Shape | Description |
|------|-------|-------------|
| `latent_cd_dae_dec.npy` | 215 × 64 | DAE-DEC optimised latent space CD |
| `latent_uc_dae_dec.npy` | 430 × 128 | DAE-DEC optimised latent space UC |

#### UK Biobank — Stability

| File | Shape | Description |
|------|-------|-------------|
| `stability_cd.npy` | (215,) | Per-patient stability scores CD (mean = 0.893, unstable = 0) |
| `stability_uc.npy` | (430,) | Per-patient stability scores UC (mean = 0.844, unstable = 0) |

#### IBDome — Preprocessed Data

| File | Shape | Description |
|------|-------|-------------|
| `ibdome_cd_preprocessed.csv` | 201 × 61 | IBDome CD protein matrix after QC |
| `ibdome_uc_preprocessed.csv` | 132 × 61 | IBDome UC protein matrix after QC |
| `ibdome_cd_final_v2.csv` | 201 × 115 | IBDome CD metadata with cluster assignments |
| `ibdome_uc_final_v2.csv` | 132 × 115 | IBDome UC metadata with cluster assignments |

#### IBDome — Cluster Labels

| File | Shape | Description |
|------|-------|-------------|
| `ibdome_labels_dec_consensus.npy` | (201,) | IBDome CD consensus labels |
| `ibdome_uc_labels_dec_consensus.npy` | (132,) | IBDome UC consensus labels |
| `ibdome_consensus_dec.npy` | 201 × 201 | IBDome CD co-occurrence matrix |
| `ibdome_uc_consensus_dec.npy` | 132 × 132 | IBDome UC co-occurrence matrix |

#### IBDome — Latent Representations

| File | Shape | Description |
|------|-------|-------------|
| `ibdome_latent_dae_dec.npy` | 201 × 16 | IBDome CD DAE-DEC latent space |
| `ibdome_uc_latent_dae_dec.npy` | 132 × 16 | IBDome UC DAE-DEC latent space |

#### IBDome — Validation Results

| File | Key Results | Description |
|------|-------------|-------------|
| `ibdome_cd_validation_v2.json` | r = 0.7346, Δ = 0.0144 | IBDome CD chemokine-vascular validation |
| `ibdome_uc_validation_v2.json` | r = 0.6886, Δ = 0.0604 | IBDome UC chemokine-vascular validation |

---

### tables/

All differential abundance,
Random Forest, pathway enrichment,
and DEC improvement results.

#### Differential Abundance

| File | Shape | Description |
|------|-------|-------------|
| `cd_limma_M3.csv` | 991 × 7 | UKB CD limma Model 3 results. Columns: logFC, AveExpr, t, P.Value, adj.P.Val, B, protein. Top: DBI logFC = −1.34, FDR = 2.12×10⁻³⁴ |
| `uc_limma_M3.csv` | 991 × 7 | UKB UC limma Model 3 results. Top: TBCB logFC = −1.74, FDR = 2.63×10⁻⁷³ |
| `ibdome_cd_de_v2.csv` | 61 × 8 | IBDome CD Welch t-test results. 58/61 significant. Top: HGF logFC = +1.33, FDR = 3.86×10⁻²¹ |
| `ibdome_uc_de_v2.csv` | 61 × 8 | IBDome UC Welch t-test results. 53/61 significant. Top: MMP-10 logFC = +1.25, FDR = 1.28×10⁻¹³ |

#### Random Forest Importance

| File | Shape | Description |
|------|-------|-------------|
| `cd_rf_importance.csv` | 991 × 2 | UKB CD Gini importance scores. Columns: protein, gini_imp. AUC = 0.983 ± 0.014 |
| `uc_rf_importance.csv` | 991 × 2 | UKB UC Gini importance scores. AUC = 0.985 ± 0.003 |
| `ibdome_cd_rf_v2.csv` | 61 × 2 | IBDome CD RF importance. AUC = 0.992 ± 0.006 |
| `ibdome_uc_rf_v2.csv` | 61 × 2 | IBDome UC RF importance. AUC = 0.998 ± 0.003 |

#### Pathway Enrichment (ORA)

| File | Pathways | Description |
|------|----------|-------------|
| `ora_cd_C1_up_Hallmark.csv` | 32 significant | UKB CD Hallmark ORA. Top: EMT FDR = 9.50×10⁻²¹ |
| `ora_uc_C1_up_Hallmark.csv` | 31 significant | UKB UC Hallmark ORA. Top: Allograft Rejection FDR = 3.36×10⁻²⁰ |
| `ibdome_cd_ora_C2_up_Hallmark.csv` | 10 significant | IBDome CD ORA. Top: TNFα/NFκB FDR = 1.32×10⁻¹⁰ |
| `ibdome_uc_ora_C2_up_Hallmark.csv` | 10 significant | IBDome UC ORA. Top: TNFα/NFκB FDR = 4.01×10⁻¹¹ |

#### DEC Improvement

| File | Shape | Description |
|------|-------|-------------|
| `dec_results_cd.csv` | 5 × 4 | CD silhouette before/after DEC for all 5 architectures. Columns: model, sil_before, sil_after, gain |
| `dec_results_uc.csv` | 5 × 4 | UC silhouette before/after DEC. Max gain: +0.580 (Standard AE) |

#### PPI Networks

| File | Description |
|------|-------------|
| `ppi_CD_C1_top50.csv` | UKB CD top 50 hub proteins with degree and betweenness centrality |
| `ppi_UC_C1_top50.csv` | UKB UC top 50 hub proteins. Top hub: SRC kinase (degree = 17, BC = 0.605) |
| `ibdome_cd_ppi_C2_top30.csv` | IBDome CD top 30 hubs. Top: IL-6, IL-10 (degree = 26) |
| `ibdome_uc_ppi_C2_top30.csv` | IBDome UC top 30 hubs |

---

### figures/

All figures generated by
`Figures.ipynb` and
`Supplementary_figures.ipynb`.

#### Main Figures

| File | Figure | Description |
|------|--------|-------------|
| `figure1_overview.png` | Figure 1 | Study design and preprocessing pipeline |
| `figure2_clustering.png` | Figure 2 | Ensemble pipeline architecture and consensus results |
| `figure3_final.png` | Figure 3 | UMAP visualisation all four cohorts |
| `figure4_final.png` | Figure 4 | Pathway enrichment and volcano plots |
| `figure5_ppi_clean.png` | Figure 5 | PPI networks hyperinflammatory subtype |

#### Supplementary Figures

| File | Figure | Description |
|------|--------|-------------|
| `figS1_ukb_cd_clinical.png` | S1 | Clinical characteristics UKB CD (n = 215) |
| `figS2_ukb_uc_clinical.png` | S2 | Clinical characteristics UKB UC (n = 430) |
| `figS3_ibdome_cd_clinical.png` | S3 | Clinical characteristics IBDome CD (n = 201) |
| `figS4_ibdome_uc_clinical.png` | S4 | Clinical characteristics IBDome UC (n = 132) |
| `figS5_consensus_matrices.png` | S5 | All four consensus co-occurrence matrices |
| `figS6_rf_bubble.png` | S6 | RF importance vs differential abundance bubble charts |
| `figS7_roc_curves.png` | S7 | ROC curves all cohorts AUC 0.983–0.998 |
| `figS9_ora_full.png` | S8 | Complete Hallmark pathway enrichment 83 pathways |
| `figS10_imputation.png` | S9 | Imputation quality assessment KS statistics |
| `figS12_stability.png` | S10 | Per-patient stability distributions |
| `figS13_dec_improvement.png` | S11 | DEC silhouette improvement all architectures |
| `figS14_loss_curves.png` | S12 | Autoencoder training and validation loss curves |
---

## Troubleshooting

ERROR: CUDA not available
Check GPU allocation:
nvidia-smi
Ensure --gres=gpu:a100:1
in SLURM script

ERROR: Module not found
Run: module purge
Then reload all modules
in the correct order

ERROR: Out of memory
Reduce batch size in notebook:
BATCH_SIZE = 16 (from 32)
Or request more RAM:
--mem=128G

ERROR: Enrichr API timeout
Check internet connection
Retry after 5 minutes
BlueBEAR nodes have outbound
internet access

ERROR: R packages not found
Ensure R module loaded:
module load R/4.3.1-foss-2023a
Reinstall in R console

ERROR: Permission denied on /rds/
Check your storage quota:
quota -s
Delete temporary files if full


---

## Citation

```bibtex
@mastersthesis{jahnavi2026prosubIBD,
    author  = {Tippanaboyana Jahnavi},
    title   = {ProSubIBD: Proteomic
               Profiling to Identify
               and Validate IBD Subtypes},
    school  = {University of Birmingham},
    year    = {2026},
    type    = {MSc Health Data Science
               Dissertation},
    supervisor = {Dr Animesh Acharjee},
    url     = {https://github.com/
               tippanaboyanajahnavi02/
               -ibd-proteomic-subtyping}
}
```

---

## Contact

**Student:** Tippanaboyana Jahnavi
**ID:** 2989054
**Programme:** MSc Health Data Science
**University:** University of Birmingham
**Supervisor:** Dr Animesh Acharjee
**Department:** Cancer and Genomic Sciences

---

## Licence

This code is released for academic
reproducibility purposes under the
MIT Licence.

UK Biobank data: access must be
obtained independently via approved
application at ukbiobank.ac.uk

IBDome data: publicly available
at ibdome.de
## Data Access

### UK Biobank
