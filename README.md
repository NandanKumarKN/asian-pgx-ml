# AI-Driven Pharmacogenomics: Predicting Drug Response Variability Across Asian Ethnic Subgroups

[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Status: In Progress](https://img.shields.io/badge/status-in%20progress-orange.svg)]()

**Author:** Nandan Kumar K N  
**Affiliation:** Independent researcher | M.Sc. Bioinformatics, University of West London  
**Contact:** nandankumarkn4@gmail.com

---

## Overview

This repository contains the full analysis pipeline for a pharmacogenomics research study investigating drug response variability across six disaggregated Asian ethnic subgroups using multi-omics integration and interpretable machine learning.

**The core question:** Does treating "Asian" as a monolithic category in pharmacogenomics guidelines mask clinically meaningful variation between South Asian (Indian, Bengali) and East Asian (Chinese, Japanese) populations?

**Target journal:** Frontiers in Pharmacology / Briefings in Bioinformatics

---

## Target pharmacogenes

| Gene | Key drugs affected | Known intra-Asian divergence |
|---|---|---|
| CYP2D6 | Codeine, tamoxifen, antidepressants | High: *10 dominant EAS; *4/*41 enriched SAS |
| CYP2C19 | Clopidogrel, PPIs, SSRIs | Very high: PM ~3% (SAS) vs ~20% (EAS) |
| CYP3A5 | Tacrolimus, immunosuppressants | Moderate: *3 null allele 15–40% |
| NUDT15 | Azathioprine, 6-mercaptopurine | Significant in EAS; rare in SAS |
| SLCO1B1 | Simvastatin (myopathy risk) | *5 allele varies across populations |

---

## Study populations (1000 Genomes Phase 3)

| Code | Population | Super-population | n |
|---|---|---|---|
| GIH | Gujarati Indian in Houston | SAS | 103 |
| ITU | Indian Telugu in the UK | SAS | 102 |
| BEB | Bengali in Bangladesh | SAS | 86 |
| CHB | Han Chinese in Beijing | EAS | 103 |
| CHS | Southern Han Chinese | EAS | 105 |
| JPT | Japanese in Tokyo | EAS | 104 |

---

## Data sources (all free, open access)

- **PharmGKB** — pharmgkb.org (variant-drug-phenotype annotations)
- **1000 Genomes Project Phase 3** — ftp.1000genomes.ebi.ac.uk (VCF files)
- **GTEx v10** — gtexportal.org (whole blood gene expression)
- **KEGG** — rest.kegg.jp (pathway annotations)

---

## Project structure

```
asian-pgx-ml/
├── data/
│   ├── raw/
│   │   ├── 1kgp/          # Downloaded VCF files (not committed — see download_data.py)
│   │   ├── pharmgkb/      # PharmGKB annotation TSVs
│   │   └── gtex/          # GTEx whole blood TPM data
│   └── processed/         # Feature matrices, allele call outputs
├── notebooks/
│   ├── 01_eda_allele_frequencies.ipynb
│   ├── 02_feature_matrix_construction.ipynb
│   ├── 03_ml_model_training.ipynb
│   ├── 04_shap_explainability.ipynb
│   └── 05_clinical_dosing_mapping.ipynb
├── src/
│   ├── download_data.py   # Fetch all raw datasets
│   ├── allele_calling.py  # PyPGx star allele pipeline
│   ├── feature_builder.py # SNP + GTEx feature matrix
│   └── train_models.py    # RF, XGBoost, Elastic Net
├── results/
│   ├── figures/           # Publication-ready plots
│   ├── tables/            # Allele frequency tables, model metrics
│   └── models/            # Saved model objects (.pkl)
├── docs/
│   ├── manuscript.docx    # Paper draft
│   └── references.bib     # BibTeX references
├── requirements.txt
├── environment.yml
├── Dockerfile
└── README.md
```

---

## Quickstart

### 1. Clone and set up environment

```bash
git clone https://github.com/YOUR_USERNAME/asian-pgx-ml.git
cd asian-pgx-ml

# Option A: conda (recommended)
conda env create -f environment.yml
conda activate pgx-ml

# Option B: pip
pip install -r requirements.txt
```

### 2. Download raw data

```bash
python src/download_data.py --all
```

### 3. Run allele calling

```bash
python src/allele_calling.py --gene CYP2C19 --populations GIH ITU BEB CHB CHS JPT
```

### 4. Build feature matrices

```bash
python src/feature_builder.py --include-gtex
```

### 5. Train models

```bash
python src/train_models.py --model all --cv 5
```

### 6. Run notebooks in order

Open Jupyter and run notebooks 01 through 05 sequentially.

---

## Analysis pipeline

```
1000 Genomes VCF + PharmGKB → PyPGx allele calling → metaboliser phenotype labels
                                                              ↓
GTEx whole blood TPM ─────────────────────────────→ feature matrix construction
                                                              ↓
                                          RF / XGBoost / Elastic Net (5-fold CV)
                                                              ↓
                                        SHAP explainability + CPIC dosing mapping
```

---

## Reproducibility

All data sources are freely available. No institutional access required.  
A Docker container is provided for full reproducibility: `docker build -t pgx-ml .`

---

## License

MIT License — see LICENSE file.

---

## Citation

If you use this code or pipeline, please cite this repository until the manuscript is published:

```
Nandan Kumar K N (2026). AI-driven pharmacogenomics: predicting drug response 
variability across Asian ethnic subgroups using multi-omics integration. 
GitHub: https://github.com/YOUR_USERNAME/asian-pgx-ml
```
