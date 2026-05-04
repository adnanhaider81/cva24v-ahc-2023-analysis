#!/usr/bin/env python3
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]


def read_fastq(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) % 4:
        raise SystemExit(f"{path} does not contain complete FASTQ records")
    records = []
    for i in range(0, len(lines), 4):
        name, seq, plus, qual = lines[i : i + 4]
        if not name.startswith("@") or plus != "+":
            raise SystemExit(f"{path} has an invalid FASTQ record near line {i + 1}")
        if len(seq) != len(qual):
            raise SystemExit(f"{path} has sequence/quality length mismatch near line {i + 1}")
        records.append((name, seq))
    return records


def main():
    cfg = yaml.safe_load((ROOT / "config" / "config.yaml").read_text(encoding="utf-8"))
    pairs = cfg.get("pairs", [])
    if not pairs:
        raise SystemExit("config/config.yaml must define at least one input pair")

    for pair in pairs:
        sample = pair["sample"]
        if sample.startswith("AHC_2023_ISB"):
            raise SystemExit("default config should use synthetic sample IDs")
        r1 = ROOT / pair["r1"]
        r2 = ROOT / pair["r2"]
        if not r1.exists() or not r2.exists():
            raise SystemExit(f"missing example FASTQ files for {sample}")
        n1 = len(read_fastq(r1))
        n2 = len(read_fastq(r2))
        if n1 != n2:
            raise SystemExit(f"paired FASTQ record count mismatch for {sample}: {n1} vs {n2}")

    contigs = ROOT / "data-example" / "contigs" / "CVA24V_DEMO_01_contigs.fasta"
    if not contigs.exists() or not contigs.read_text(encoding="utf-8").startswith(">"):
        raise SystemExit("synthetic contig FASTA is missing or invalid")

    print("Synthetic example inputs validated")


if __name__ == "__main__":
    main()
