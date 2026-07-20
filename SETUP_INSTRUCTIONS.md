# Setup Instructions — Run these commands once

## Step 1: Install Git (if not already installed)
Download from: https://git-scm.com/downloads

## Step 2: Create your GitHub repo
1. Go to github.com → New repository
2. Name it: asian-pgx-ml
3. Set to PUBLIC (important for PhD applications — professors will look at it)
4. Do NOT initialise with README (you already have one)
5. Copy the repo URL

## Step 3: Initialise Git locally
Open terminal in the project folder and run:

```bash
git init
git add .
git commit -m "Initial project structure: pharmacogenomics ML pipeline"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/asian-pgx-ml.git
git push -u origin main
```

## Step 4: Set up the Python environment
```bash
# Install Miniconda first if you don't have it:
# https://docs.conda.io/en/latest/miniconda.html

conda env create -f environment.yml
conda activate pgx-ml

# Verify key packages installed:
python -c "import pandas, numpy, sklearn, xgboost, shap; print('All packages OK')"
```

## Step 5: Download the data
```bash
python src/download_data.py --pharmgkb   # Start with PharmGKB (smallest)
python src/download_data.py --1kgp       # Generates shell script for VCF download
python src/download_data.py --gtex       # ~1.5 GB — run overnight
```

## Step 6: Open Jupyter and start notebook 01
```bash
jupyter notebook notebooks/01_eda_allele_frequencies.ipynb
```

---
You're now set up. Focus on Month 1: read 10 papers and explore the PharmGKB data.
