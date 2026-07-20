"""
download_data.py
----------------
Downloads all raw data required for the asian-pgx-ml study:
  1. PharmGKB annotation tables (VIP genes, clinical annotations)
  2. 1000 Genomes Project Phase 3 VCF files (target pharmacogene regions)
  3. GTEx v10 whole blood gene expression data

Usage:
    python src/download_data.py --all
    python src/download_data.py --pharmgkb
    python src/download_data.py --1kgp --genes CYP2D6 CYP2C19
    python src/download_data.py --gtex
"""

import argparse
import os
import sys
import requests
from pathlib import Path
from tqdm import tqdm

# ── Project root (one level up from src/) ──────────────────────────────────
ROOT = Path(__file__).resolve().parent.parent
DATA_RAW = ROOT / "data" / "raw"


# ── Target pharmacogenes and their GRCh38 coordinates ──────────────────────
GENE_COORDS = {
    "CYP2D6":  ("chr22", 42_126_499, 42_130_810),
    "CYP2C19": ("chr10", 94_762_681, 94_855_547),
    "CYP3A5":  ("chr7",  99_245_157, 99_277_590),
    "NUDT15":  ("chr13", 48_037_168, 48_041_198),
    "SLCO1B1": ("chr12", 21_282_128, 21_395_730),
}

# ── 1KGP populations to include ─────────────────────────────────────────────
POPULATIONS = ["GIH", "ITU", "BEB", "CHB", "CHS", "JPT"]

# ── PharmGKB download URLs ───────────────────────────────────────────────────
PHARMGKB_URLS = {
    "relationships":         "https://api.pharmgkb.org/v1/download/file/data/relationships.zip",
    "clinical_annotations":  "https://api.pharmgkb.org/v1/download/file/data/clinicalAnnotations.zip",
    "var_drug_ann":          "https://api.pharmgkb.org/v1/download/file/data/variantAnnotations.zip",
    "genes":                 "https://api.pharmgkb.org/v1/download/file/data/genes.zip",
}

# ── GTEx v10 whole blood TPM ─────────────────────────────────────────────────
GTEX_URL = (
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/"
    "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz"
)


def download_file(url: str, dest: Path, label: str = "") -> bool:
    """Stream-download a file with a tqdm progress bar."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        print(f"  [skip] {dest.name} already exists")
        return True
    try:
        r = requests.get(url, stream=True, timeout=60)
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        desc = label or dest.name
        with open(dest, "wb") as f, tqdm(
            total=total, unit="B", unit_scale=True, desc=desc[:40]
        ) as bar:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
                bar.update(len(chunk))
        print(f"  [ok] saved to {dest}")
        return True
    except Exception as e:
        print(f"  [error] {label}: {e}")
        if dest.exists():
            dest.unlink()
        return False


def download_pharmgkb():
    """Download PharmGKB annotation files."""
    print("\n── PharmGKB ─────────────────────────────────")
    out = DATA_RAW / "pharmgkb"
    out.mkdir(parents=True, exist_ok=True)
    for name, url in PHARMGKB_URLS.items():
        download_file(url, out / f"{name}.zip", label=name)
    print("  Done. Unzip files manually or add --unzip flag.")


def download_1kgp(genes: list = None):
    """
    Download 1KGP Phase 3 VCF slices for target pharmacogene regions.
    Uses bcftools view to pull only the relevant genomic windows.

    NOTE: This function constructs the commands. Run them in your shell
    after installing bcftools (conda install -c bioconda bcftools).
    """
    genes_to_use = genes or list(GENE_COORDS.keys())
    out = DATA_RAW / "1kgp"
    out.mkdir(parents=True, exist_ok=True)

    base_url = (
        "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
        "1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/"
    )

    print("\n── 1000 Genomes Project Phase 3 ─────────────────────────────")
    print("  Generating bcftools commands for pharmacogene region extraction...\n")

    commands = []
    for gene in genes_to_use:
        if gene not in GENE_COORDS:
            print(f"  [warn] Unknown gene: {gene}")
            continue
        chrom, start, end = GENE_COORDS[gene]
        # Add ±5kb flanking for star allele calling accuracy
        region = f"{chrom}:{max(1, start-5000)}-{end+5000}"
        vcf_name = chrom.replace("chr", "")
        remote_vcf = (
            f"{base_url}ALL.{vcf_name}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
        )
        out_file = out / f"{gene}_GRCh38.vcf.gz"
        pop_filter = "|".join(POPULATIONS)

        cmd = (
            f"bcftools view "
            f'  --regions "{region}" '
            f'  --samples-file <(bcftools query -l {remote_vcf} | grep -E "{pop_filter}") '
            f'  --output-type z '
            f'  --output {out_file} '
            f'  "{remote_vcf}" && '
            f"bcftools index --tbi {out_file}"
        )
        commands.append((gene, cmd))
        print(f"  [{gene}] region: {region}")
        print(f"  Command:\n    {cmd}\n")

    # Save commands to a shell script
    script_path = out / "run_1kgp_download.sh"
    with open(script_path, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("# Auto-generated by download_data.py\n")
        f.write("# Run this script after installing bcftools\n\n")
        for gene, cmd in commands:
            f.write(f"echo 'Downloading {gene}...'\n")
            f.write(cmd + "\n\n")
    print(f"  Shell script saved to: {script_path}")
    print("  Run it with: bash data/raw/1kgp/run_1kgp_download.sh")


def download_gtex():
    """
    Download GTEx v10 whole blood gene expression (TPM matrix).
    WARNING: This file is ~1.5 GB. Make sure you have disk space.
    """
    print("\n── GTEx v10 whole blood ─────────────────────────────────────")
    out = DATA_RAW / "gtex"
    out.mkdir(parents=True, exist_ok=True)
    dest = out / "GTEx_v10_whole_blood_tpm.gct.gz"
    print("  NOTE: GTEx TPM file is ~1.5 GB. This will take several minutes.")
    print(f"  Downloading to: {dest}")
    download_file(GTEX_URL, dest, label="GTEx v10 whole blood TPM")
    print("\n  After download, extract pharmacogene rows with:")
    genes_str = "|".join(GENE_COORDS.keys())
    print(f'    zcat {dest} | grep -E "^#{genes_str}" > data/raw/gtex/pharmacogenes_tpm.tsv')


def main():
    parser = argparse.ArgumentParser(
        description="Download raw data for asian-pgx-ml study"
    )
    parser.add_argument("--all", action="store_true", help="Download all datasets")
    parser.add_argument("--pharmgkb", action="store_true", help="Download PharmGKB files")
    parser.add_argument("--1kgp", dest="kgp", action="store_true", help="Generate 1KGP download commands")
    parser.add_argument("--gtex", action="store_true", help="Download GTEx data")
    parser.add_argument("--genes", nargs="+", default=None,
                        help="Genes to include for 1KGP (default: all 5 pharmacogenes)")
    args = parser.parse_args()

    if not any([args.all, args.pharmgkb, args.kgp, args.gtex]):
        parser.print_help()
        sys.exit(0)

    print(f"\nProject root: {ROOT}")
    print(f"Output dir:   {DATA_RAW}\n")

    if args.all or args.pharmgkb:
        download_pharmgkb()
    if args.all or args.kgp:
        download_1kgp(genes=args.genes)
    if args.all or args.gtex:
        download_gtex()

    print("\n── All done ──────────────────────────────────────────────────")
    print("Next step: python src/allele_calling.py --gene CYP2C19")


if __name__ == "__main__":
    main()
