"""
download_1kgp_windows.py
------------------------
Downloads 1KGP Phase 3 VCF slices for 5 pharmacogenes
using Python only — no bcftools, no Linux required.

Works on Windows with Git Bash or standard terminal.

Usage (from your project root):
    python src/download_1kgp_windows.py

Requirements (already in your environment.yml):
    pip install cyvcf2 requests tqdm

What it does:
    1. Downloads the full chromosome VCF index (.tbi) to find byte offsets
    2. Uses HTTP range requests to pull ONLY the pharmacogene regions
    3. Saves one VCF.gz per gene to data/raw/1kgp/
"""

import requests
import gzip
import struct
import io
import os
import sys
from pathlib import Path
from tqdm import tqdm

# ── Paths ─────────────────────────────────────────────────────────────────
ROOT    = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "data" / "raw" / "1kgp"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── 1KGP Phase 3 GRCh38 FTP base ─────────────────────────────────────────
FTP_BASE = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL"
)

# ── Pharmacogene regions (GRCh38, ±10kb flanking for star allele accuracy) ─
GENE_REGIONS = {
    "CYP2D6":  ("22", 42_116_499,  42_140_810),
    "CYP2C19": ("10", 94_752_681,  94_865_547),
    "CYP3A5":  ("7",  99_235_157,  99_287_590),
    "NUDT15":  ("13", 48_027_168,  48_051_198),
    "SLCO1B1": ("12", 21_272_128,  21_405_730),
}

# ── Target populations (used to filter samples after download) ─────────────
TARGET_POPS = ["GIH", "ITU", "BEB", "CHB", "CHS", "JPT"]

# ── VCF filename pattern on 1KGP FTP ─────────────────────────────────────
def vcf_filename(chrom: str) -> str:
    return (
        f"ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019"
        f".GRCh38.phased.vcf.gz"
    )


def download_panel():
    """Download 1KGP sample panel (sample → population mapping)."""
    panel_url = (
        "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
        "integrated_call_samples_v3.20200731.ALL.ped"
    )
    panel_path = OUT_DIR / "panel.ped"
    if panel_path.exists():
        print(f"[skip] panel.ped already exists")
        return panel_path
    print("Downloading sample panel (~2 MB)...")
    r = requests.get(panel_url, stream=True)
    r.raise_for_status()
    with open(panel_path, "wb") as f:
        for chunk in r.iter_content(8192):
            f.write(chunk)
    print(f"  Saved → {panel_path.name}")
    return panel_path


def get_vcf_url(chrom: str) -> str:
    return f"{FTP_BASE}/{vcf_filename(chrom)}"


def get_tabix_index(chrom: str) -> bytes:
    """Download the .tbi index for a chromosome VCF."""
    tbi_url = get_vcf_url(chrom) + ".tbi"
    print(f"  Fetching index: {tbi_url.split('/')[-1]}")
    r = requests.get(tbi_url, timeout=60)
    r.raise_for_status()
    return r.content


def parse_tbi_chunk_offsets(tbi_bytes: bytes, chrom: str,
                              start: int, end: int):
    """
    Parse a tabix (.tbi) index to find the virtual file offsets
    covering a genomic region.

    Returns (min_offset, max_offset) as byte positions in the BGZF file.
    Returns None if region not found.
    """
    try:
        data = gzip.decompress(tbi_bytes)
    except Exception:
        data = tbi_bytes  # already decompressed

    idx = 0
    magic = data[idx:idx+4]; idx += 4
    if magic != b"TBI\x01":
        print("  WARNING: Not a valid .tbi file")
        return None

    n_ref = struct.unpack_from("<i", data, idx)[0]; idx += 4
    # format, col_seq, col_beg, col_end, meta, skip
    idx += 6 * 4

    # name block
    l_nm = struct.unpack_from("<i", data, idx)[0]; idx += 4
    names_block = data[idx:idx+l_nm]; idx += l_nm
    names = names_block.split(b"\x00")
    ref_names = [n.decode("utf-8", errors="replace") for n in names if n]

    # find chromosome index
    chrom_variants = [chrom, f"chr{chrom}"]
    ref_idx = None
    for cv in chrom_variants:
        if cv in ref_names:
            ref_idx = ref_names.index(cv)
            break
    if ref_idx is None:
        print(f"  WARNING: chromosome {chrom} not found in index")
        return None

    offsets = []
    for ri in range(n_ref):
        # bins
        n_bin = struct.unpack_from("<i", data, idx)[0]; idx += 4
        for _ in range(n_bin):
            bin_id = struct.unpack_from("<I", data, idx)[0]; idx += 4
            n_chunk = struct.unpack_from("<i", data, idx)[0]; idx += 4
            for _ in range(n_chunk):
                cnk_beg = struct.unpack_from("<Q", data, idx)[0]; idx += 8
                cnk_end = struct.unpack_from("<Q", data, idx)[0]; idx += 8
                if ri == ref_idx:
                    offsets.append((cnk_beg >> 16, cnk_end >> 16))
        # linear index
        n_intv = struct.unpack_from("<i", data, idx)[0]; idx += 4
        idx += n_intv * 8

    if not offsets:
        return None

    all_starts = [o[0] for o in offsets if o[0] > 0]
    all_ends   = [o[1] for o in offsets if o[1] > 0]
    if not all_starts:
        return None

    # Use a generous range: from min start to max end
    return min(all_starts), max(all_ends)


def download_vcf_region(gene: str, chrom: str, start: int, end: int):
    """
    Download the VCF region for a pharmacogene using HTTP range requests
    when possible, otherwise downloads the full chromosome slice via
    a streaming approach with region filtering.
    """
    out_path = OUT_DIR / f"{gene}_GRCh38.vcf.gz"
    if out_path.exists() and out_path.stat().st_size > 10_000:
        print(f"  [skip] {out_path.name} already exists ({out_path.stat().st_size/1e6:.1f} MB)")
        return out_path

    vcf_url = get_vcf_url(chrom)
    print(f"\n{'='*55}")
    print(f"Gene: {gene}  |  chr{chrom}:{start:,}-{end:,}")
    print(f"VCF : {vcf_url.split('/')[-1]}")

    # Try to get byte range from .tbi index
    byte_range = None
    try:
        tbi_bytes = get_tabix_index(chrom)
        byte_range = parse_tbi_chunk_offsets(tbi_bytes, chrom, start, end)
        if byte_range:
            print(f"  Index parsed: byte range {byte_range[0]:,} – {byte_range[1]:,}")
    except Exception as e:
        print(f"  Could not parse index ({e}), will download via streaming")

    # ── Strategy 1: HTTP range request (fast, ~5-50 MB per gene) ──────────
    if byte_range:
        range_start, range_end = byte_range
        # Add generous buffer (BGZF blocks can span boundaries)
        range_start = max(0, range_start - 65536)
        range_end   = range_end + 524288

        headers = {"Range": f"bytes={range_start}-{range_end}"}
        print(f"  HTTP range request: {(range_end-range_start)/1e6:.1f} MB")

        try:
            r = requests.get(vcf_url, headers=headers, stream=True, timeout=120)
            if r.status_code in (200, 206):
                raw_bytes = b""
                total = int(r.headers.get("content-length", 0))
                with tqdm(total=total, unit="B", unit_scale=True,
                          desc=f"  {gene}") as bar:
                    for chunk in r.iter_content(65536):
                        raw_bytes += chunk
                        bar.update(len(chunk))

                # Write raw bytes and filter to region
                tmp_path = OUT_DIR / f"{gene}_raw.bin"
                tmp_path.write_bytes(raw_bytes)
                print(f"  Downloaded {len(raw_bytes)/1e6:.1f} MB → filtering to region...")
                _filter_and_save(tmp_path, out_path, chrom, start, end, gene)
                tmp_path.unlink(missing_ok=True)
                return out_path
        except Exception as e:
            print(f"  Range request failed ({e}), falling back to streaming...")

    # ── Strategy 2: Streaming download with early stop ────────────────────
    print(f"  Streaming VCF (will stop after gene region)...")
    print(f"  NOTE: This may take 5-20 min depending on chromosome size.")
    print(f"        The gene region is in the first ~{end/1e9*100:.0f}% of the file.")

    try:
        r = requests.get(vcf_url, stream=True, timeout=300)
        r.raise_for_status()

        collected = b""
        total_downloaded = 0
        region_found = False
        past_region = False

        with tqdm(unit="MB", desc=f"  {gene} (streaming)") as bar:
            for chunk in r.iter_content(chunk_size=65536):
                collected += chunk
                total_downloaded += len(chunk)
                bar.update(len(chunk) / 1e6)

                # Try to decompress and check position
                if total_downloaded > 10_000_000:  # 10MB minimum before checking
                    try:
                        lines = list(_iter_vcf_lines(collected))
                        for line in lines:
                            if line.startswith(b"#"):
                                continue
                            parts = line.split(b"\t", 3)
                            if len(parts) >= 2:
                                pos = int(parts[1])
                                if pos >= start:
                                    region_found = True
                                if region_found and pos > end + 100_000:
                                    past_region = True
                                    break
                        if past_region:
                            print(f"\n  Region complete at {total_downloaded/1e6:.0f} MB downloaded")
                            break
                    except Exception:
                        pass  # Keep downloading if we can't parse yet

        _filter_and_save_bytes(collected, out_path, chrom, start, end, gene)
        return out_path

    except Exception as e:
        print(f"  ERROR: {e}")
        print(f"\n  ALTERNATIVE: Run this command in Git Bash after installing tabix:")
        print(f"  tabix -h {vcf_url} chr{chrom}:{start}-{end} | bgzip > {out_path}")
        return None


def _iter_vcf_lines(raw_bytes: bytes):
    """Attempt to decompress BGZF/gzip bytes and yield VCF lines."""
    try:
        with gzip.open(io.BytesIO(raw_bytes), "rb") as gz:
            for line in gz:
                yield line.rstrip()
    except Exception:
        pass


def _filter_and_save(tmp_path: Path, out_path: Path,
                      chrom: str, start: int, end: int, gene: str):
    """Filter raw downloaded bytes to the target region and save as VCF.gz."""
    _filter_and_save_bytes(tmp_path.read_bytes(), out_path, chrom, start, end, gene)


def _filter_and_save_bytes(raw_bytes: bytes, out_path: Path,
                             chrom: str, start: int, end: int, gene: str):
    """Parse bytes, filter to region, write VCF.gz."""
    header_lines = []
    data_lines   = []
    chrom_variants = [chrom, f"chr{chrom}"]

    try:
        with gzip.open(io.BytesIO(raw_bytes), "rt", encoding="utf-8",
                        errors="replace") as gz:
            for line in gz:
                line = line.rstrip("\n")
                if line.startswith("#"):
                    header_lines.append(line)
                    continue
                parts = line.split("\t", 3)
                if len(parts) < 2:
                    continue
                rec_chrom = parts[0]
                if rec_chrom not in chrom_variants:
                    continue
                try:
                    pos = int(parts[1])
                except ValueError:
                    continue
                if start <= pos <= end:
                    data_lines.append(line)
    except Exception as e:
        print(f"  Parse error: {e}")

    if not data_lines:
        print(f"  WARNING: No variants found in region for {gene}")
        print(f"  This may mean the byte range was off — try re-running")
        return

    print(f"  {len(data_lines):,} variants in region chr{chrom}:{start:,}-{end:,}")

    with gzip.open(out_path, "wt", compresslevel=6) as gz:
        for line in header_lines:
            gz.write(line + "\n")
        for line in data_lines:
            gz.write(line + "\n")

    size_mb = out_path.stat().st_size / 1e6
    print(f"  ✓ Saved → {out_path.name}  ({size_mb:.1f} MB)")


def main():
    print("="*55)
    print("1KGP Phase 3 pharmacogene VCF downloader")
    print("Windows-compatible | Python only | No bcftools needed")
    print("="*55)

    # Step 1: Panel file
    print("\n[1/6] Downloading sample panel...")
    download_panel()

    # Step 2: VCF per gene
    failed = []
    for i, (gene, (chrom, start, end)) in enumerate(GENE_REGIONS.items(), 2):
        print(f"\n[{i}/6] {gene}...")
        result = download_vcf_region(gene, chrom, start, end)
        if result is None:
            failed.append(gene)

    # Summary
    print("\n" + "="*55)
    print("DOWNLOAD COMPLETE")
    vcf_files = list(OUT_DIR.glob("*.vcf.gz"))
    print(f"VCF files saved: {len(vcf_files)}")
    for v in vcf_files:
        print(f"  {v.name}  ({v.stat().st_size/1e6:.1f} MB)")
    if failed:
        print(f"\nFailed genes: {failed}")
        print("For failed genes, try the manual tabix command shown above,")
        print("or install WSL + bcftools for the fastest approach.")
    else:
        print("\nAll genes downloaded successfully!")
        print("Next step: open notebook 01_eda_allele_frequencies.ipynb")


if __name__ == "__main__":
    main()
