# Genomic characterization of the Coxsackievirus A24 variant in the Acute Hemorrhagic Conjunctivitis outbreak 2023 in Islamabad, Pakistan through metagenomic next generation sequencing

Reproducible code and workflow that mirror the analysis in the manuscript.  
Primary paper DOI: [10.1016/j.jviromet.2025.115213](https://doi.org/10.1016/j.jviromet.2025.115213)

## Program summary
End to end analysis of CV-A24v using metagenomic NGS for outbreak investigation. Steps match the study design and are fully scripted so a reviewer can reproduce the workflow.

1) Inputs
   - Paired-end FASTQ from conjunctival swabs sequenced on Illumina MiSeq 2x150.

2) Quality control and trimming
   - FastQC for initial QC.
   - Trimmomatic for adapter and quality trimming.
   - Picard MarkDuplicates to remove PCR duplicates.

3) De novo assembly and contig validation
   - Assemble with SPAdes.
   - Contig QC: length, percent Ns, GC, N50; expect principal CV-A24v contig near 7.4 kb.
   - BLAST remote against NCBI nt to identify closest matches.
   - Download selected GenBank sequences for context.

4) Reference-based mapping and masked consensus
   - Map to recent 2023 CV-A24v genome or to the prototype D90457.1.
   - Sort, mark duplicates, compute depth, mask sites below a coverage threshold to N.
   - Call variants with bcftools, build masked consensus with bcftools consensus.

5) Phylogeny for whole genome and VP1
   - Whole genome alignment with MAFFT.
   - VP1 extraction, trimming to the region used in the manuscript.
   - IQ-TREE with 1000 ultrafast bootstraps and model control: GTR+G+I for whole genome and K2+I for VP1 to match the study, or MFP for auto model selection.
   - Trees saved in Newick and log files.

6) Mutation and comparative analysis
   - Open reading frame extraction and AA comparison vs the prototype strain D90457.1.
   - VP1 and VP3 comparison against 2005 Pakistan sequences AB365074 to AB365078.
   - Optional regional recombination scan by slicing 5'UTR, 2A to 2C, and 3A to 3D regions and BLASTing each separately.

7) Optional structural modeling preparation
   - Scripts generate FASTA segments ready for Robetta and guidance for UCSF Chimera to reproduce the VP1 structural comparison.

8) Outputs
   - Per sample consensus: results/consensus/<sample>.fa
   - Combined consensus: results/consensus/all_consensus.fasta
   - Whole genome alignment: results/aln/wg_alignment.fasta
   - VP1 alignment: results/aln/vp1_alignment.fasta
   - Phylogeny: results/iqtree/wg.treefile and results/iqtree/vp1.treefile
   - QC, BLAST tables, recombination segment hits, and mutation reports under results/

9) Repro and compliance
   - CITATION.cff, MIT LICENSE, CI check for quick verification.
   - Never commit restricted or clinical data. Set NCBI_EMAIL once for Entrez-based steps.

## Requirements
- Python 3.11 or newer
- Option A: pip and virtualenv
- Option B: conda or mamba
- Snakemake if you want to run the full pipeline
- Optional: Docker

### NCBI usage note
Set a contact email once per shell for E-utilities. Optional API key improves rate limits.
```bash
export NCBI_EMAIL="you@example.com"
export NCBI_API_KEY="xxxxxxxxxxxxxxxxxxxxxxxxxxxx"   # optional
```

## Quick verification
```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
python -m pip install -r env/requirements.txt
python analysis/scripts/example_qc_plot.py --in data-example/example_counts.tsv --out results-example/example_plot.png
```

## One-command end to end run
```bash
export NCBI_EMAIL="you@example.com"
conda env create -f env/environment.yml
conda activate cva24-env
snakemake -s workflow/Snakefile -c 4 --printshellcmds
```

## Configuration
Edit config/config.yaml. Example:
```yaml
pairs:
  - sample: AHC_001
    r1: data-private/AHC_001_R1.fastq.gz
    r2: data-private/AHC_001_R2.fastq.gz

ref_accessions:
  prototype: D90457.1               # CV-A24v prototype, Singapore 1970
  recent_2023: OR361387.1           # 2023 Zhongshan China complete genome

context_accessions:                 # Used to recreate study trees
  - OR361388.1
  - OR803779.1
  - OR361389.1
  - OR361390.1
  - KR478685.1
  - KR399988.1
  - AB365074.1
  - AB365075.1
  - AB365076.1
  - AB365077.1
  - AB365078.1
  - OR633288.1
  - PP028462.1
  - PP028461.1

params:
  threads: 4
  trim_adapters: env/TruSeq3-PE.fa
  min_depth_consensus: 10
  min_qual: 20
  iqtree_model_wg: GTR+G+I
  iqtree_model_vp1: K2+I
  bootstrap: 1000
```

Override models at runtime if you want auto selection:
```bash
snakemake -s workflow/Snakefile -c 4 --config params:iqtree_model_wg=MFP params:iqtree_model_vp1=MFP
```

## Structural modeling notes
- Robetta server for homology modeling of VP1.
- UCSF Chimera for alignment and surface or volume analysis.
- The repo includes a script that exports VP1 sequences for upload and a short how-to in docs/STRUCTURE.md.

## Methods and tool references
FastQC v0.11.9, Trimmomatic v0.39, Picard MarkDuplicates, SPAdes v3.15.5, BWA, MAFFT, IQ-TREE v2 with 1000 bootstraps, MEGA model selection references K2+I and GTR+G+I, BLAST+, Entrez Direct.
