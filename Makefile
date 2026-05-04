SHELL := /bin/bash

ENV_NAME := cva24v-env
THREADS ?= 4
PYTHON ?= python3

.PHONY: env run dry smoke test clean help

help:
	@echo "make env    - create conda env"
	@echo "make run    - run full Snakemake workflow"
	@echo "make dry    - dry run to show planned steps"
	@echo "make smoke  - validate bundled synthetic example inputs"
	@echo "make test   - alias for make smoke"
	@echo "make clean  - remove work and generated results"

env:
	conda env create -f env/environment.yml || echo "Env may already exist"
	@echo "Activate with: conda activate $(ENV_NAME)"

run:
	@if [ -z "$$NCBI_EMAIL" ]; then echo "Set NCBI_EMAIL before running"; exit 1; fi
	snakemake -s workflow/Snakefile -c $(THREADS) --printshellcmds

dry:
	snakemake -s workflow/Snakefile -n -c $(THREADS)

smoke:
	$(PYTHON) -m py_compile analysis/scripts/*.py
	$(PYTHON) tests/validate_example_inputs.py
	$(PYTHON) analysis/scripts/contig_qc.py --in data-example/contigs/CVA24V_DEMO_01_contigs.fasta --min_len 80 --out_tsv results-example/contig_qc.tsv
	$(PYTHON) analysis/scripts/example_qc_plot.py --in data-example/example_counts.tsv --out results-example/example_plot.png
	@echo "Wrote results-example/example_plot.png"

test: smoke

clean:
	rm -rf work results logs refs .snakemake analysis/scripts/__pycache__ tests/__pycache__ results-example/contig_qc.tsv results-example/contig_qc.summary.txt
