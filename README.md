# Genomic characterization of the Coxsackievirus A24 variant in the Acute Hemorrhagic Conjunctivitis outbreak 2023 in Islamabad, Pakistan through metagenomic next generation sequencing

Reproducible pathogen discovery and targeted analysis that mirror the Journal of Virological Methods paper. DOI: 10.1016/j.jviromet.2025.115213

## Program summary
One pipeline under Snakemake. Discovery first, then targeted CV-A24v analysis.

1) Inputs
   - Paired end FASTQ from metagenomic RNA libraries sequenced on Illumina MiSeq 2x150.

2) Quality control and trimming
   - FastQC for initial QC.
   - Trimmomatic with sliding window 4:30, minimum length 50, and leading and trailing quality 3.
   - Optional duplicate flagging after mapping with Picard MarkDuplicates.

3) Pathogen discovery
   - De novo contigs with SPAdes.
   - Kraken2 classification on reads and contigs using a local database that you supply.
   - Kaiju protein level classification of contigs using a local kaiju database.
   - Cross check dominant taxa with NCBI BLASTN in remote mode against nt, or against a local nt database.
   - Write a combined taxonomy summary and flag suggested targets for downstream analysis.

4) Targeted CV-A24v analysis
   - Select reference for mapping from BLAST top hits or fallback to a configured accession.
   - BWA MEM mapping, SAMtools sort and index.
   - Depth mask at minimum coverage threshold default 10 and variant calling with bcftools. Write masked consensus per sample.
   - Context set: fetch additional CV-A24v accessions and prototype D90457.1 from GenBank for comparison.
   - Whole genome phylogeny with MAFFT and IQ-TREE. Default model GTR+G+I and 1000 ultrafast bootstraps with optional ModelFinder.
   - VP1 region analysis: fetch annotation for the selected reference, extract VP1 coordinates, cut VP1 from consensuses, align, and infer a VP1 tree with a configurable model default K2+I.
   - Mutation summary: translate ORFs, compare to prototype D90457.1 and 2005 Pakistan strains AB365074 to AB365078 and write a table.

5) Outputs
   - `results/discovery/kraken2/<sample>.reads.report.txt` and `.contigs.report.txt`
   - `results/discovery/kaiju/<sample>.kaiju.report.tsv`
   - `results/discovery/<sample>_taxonomy_summary.tsv` and `<sample>_selected_targets.yaml`
   - `results/consensus/<sample>.fa` and combined `results/consensus/all_consensus.fasta`
   - `results/aln/wg_alignment.fasta`, `results/iqtree/wg.treefile`
   - `results/vp1/vp1_alignment.fasta`, `results/iqtree/vp1.treefile`
   - `results/mutations/aa_changes.tsv` and `results/mutations/vp1_snp_table.tsv`

## Requirements
- Python 3.11 or newer
- Option A: pip and virtualenv
- Option B: conda or mamba
- System tools installed through conda: fastqc, trimmomatic, spades, blast, bwa, samtools, bcftools, mafft, iqtree, entrez-direct, seqkit, snakemake, biopython, pyyaml

### NCBI usage note
Set a contact email once per shell for E-utilities. Optional API key improves rate limits.
```bash
export NCBI_EMAIL="you@example.com"
export NCBI_API_KEY="xxxxxxxxxxxxxxxxxxxxxxxxxxxx"   # optional
```

### Kraken2 and Kaiju databases
This repo does not ship databases. You must download and point the config at your local copies.

Authoritative sources
- Kraken2 homepage and docs: https://ccb.jhu.edu/software/kraken2/
- Kraken2 GitHub: https://github.com/DerrickWood/kraken2
- Ben Langmead AWS-hosted Kraken2 plus Bracken indexes: https://benlangmead.github.io/aws-indexes/k2
- Kaiju homepage: https://kaiju.binf.ku.dk/
- Kaiju prebuilt indexes: https://bioinformatics-centre.github.io/kaiju/downloads.html
- NCBI Taxonomy FTP for names.dmp and nodes.dmp: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/

Quick start examples
```bash
# Kraken2 database
mkdir -p /data/db/kraken2 && cd /data/db/kraken2
# choose a tarball from the Langmead page and extract here
tar -xzf k2_standard_202507.tar.gz

# Kaiju database
mkdir -p /data/db/kaiju && cd /data/db/kaiju
# choose a prebuilt archive and extract here
tar -xzf kaiju_db_nr_euk_2023-05-16.tgz
```

Configure paths in `config/config.yaml`:
```yaml
kraken2:
  enable: true
  db: /data/db/kraken2/k2_standard_202507

kaiju:
  enable: true
  db_fmi: /data/db/kaiju/kaiju_db_nr_euk.fmi
  nodes: /data/db/kaiju/nodes.dmp
  names: /data/db/kaiju/names.dmp
```

## One command run
```bash
export NCBI_EMAIL="you@example.com"
conda env create -f env/environment.yml
conda activate cva24v-env
snakemake -s workflow/Snakefile -c 4 --printshellcmds
```

## Configuration
Edit `config/config.yaml`. Minimal example:
```yaml
pairs:
  - sample: AHC_2023_ISB_01
    r1: data-private/AHC_2023_ISB_01_R1.fastq.gz
    r2: data-private/AHC_2023_ISB_01_R2.fastq.gz

discovery:
  min_len_contig: 300
  blast_remote: true

reference:
  fallback_acc: D90457.1    # prototype Singapore 1970

context_accessions:
  - D90457.1
  - AB365074.1
  - AB365075.1
  - AB365076.1
  - AB365077.1
  - AB365078.1

phylogeny:
  wg_model: GTR+G+I
  vp1_model: K2P+I
  bootstrap: 1000
  use_model_finder: false
  min_depth_consensus: 10

kraken2:
  enable: true
  db: /data/db/kraken2/k2_standard_202507

kaiju:
  enable: true
  db_fmi: /data/db/kaiju/kaiju_db_nr_euk.fmi
  nodes: /data/db/kaiju/nodes.dmp
  names: /data/db/kaiju/names.dmp
```

## How to cite
- Paper: Haider SA, Jamal Z, Ammar M, Hakim R, Afrough B, Kreku A, Inamdar L, Salman M, Umair M. Genomic characterization of the Coxsackievirus A24 variant in the Acute Hemorrhagic Conjunctivitis outbreak 2023 in Islamabad, Pakistan through metagenomic next generation sequencing. Journal of Virological Methods. 2025. https://doi.org/10.1016/j.jviromet.2025.115213
- Software: Haider SA. CV-A24v AHC 2023 Pakistan analysis. Version 2.0.0. GitHub repository.

## References
- Andrews S. 2010. FastQC. Babraham Bioinformatics.
- Bolger AM, Lohse M, Usadel B. 2014. Trimmomatic. Bioinformatics 30:2114-2120.
- Bankevich A, Nurk S, Antipov D, et al. 2012. SPAdes. J Comput Biol 19:455-477.
- Camacho C, et al. 2009. BLAST+. BMC Bioinformatics 10:421.
- Li H, 2013. BWA-MEM. arXiv:1303.3997.
- Li H, et al. 2009. SAMtools. Bioinformatics 25:2078-2079.
- Danecek P, et al. 2021. BCFtools. GigaScience 10:giab008.
- Katoh K, Standley DM. 2013. MAFFT. Mol Biol Evol 30:772-780.
- Minh BQ, et al. 2020. IQ-TREE 2. Mol Biol Evol 37:1530-1534.
- Kans J. Entrez E-utilities Help. NCBI.
- Köster J, Rahmann S. 2012. Snakemake. Bioinformatics 28:2520-2522.
- Wood DE, Lu J, Langmead B. 2019. Kraken 2. Genome Biology 20:257.
- Menzel P, Ng KL, Krogh A. 2016. Kaiju. Nat Commun 7:11257.
