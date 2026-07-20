FROM continuumio/miniconda3:24.1.2-0

LABEL maintainer="nandankumarkn4@gmail.com"
LABEL description="Reproducible environment for asian-pgx-ml pharmacogenomics study"

# Install system dependencies (bcftools, bedtools, samtools)
RUN apt-get update && apt-get install -y \
    bcftools \
    bedtools \
    samtools \
    wget \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /workspace

# Copy environment file and install
COPY environment.yml .
RUN conda env create -f environment.yml && conda clean -afy

# Make conda env default
SHELL ["conda", "run", "-n", "pgx-ml", "/bin/bash", "-c"]

# Copy project files
COPY . .

# Expose Jupyter port
EXPOSE 8888

# Default: start Jupyter
CMD ["conda", "run", "--no-capture-output", "-n", "pgx-ml", \
     "jupyter", "notebook", "--ip=0.0.0.0", "--port=8888", \
     "--no-browser", "--allow-root", "--notebook-dir=/workspace"]
