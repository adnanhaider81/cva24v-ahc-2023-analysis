# Contributing

Thank you for your interest in improving this analysis.

## Ground rules
- Do not commit patient data or restricted datasets.
- Keep parameters in `config/config.yaml` so runs are reproducible.
- Match tool versions in `env/environment.yml` where possible.

## Workflow
1. Fork the repo and create a feature branch.
2. Make focused changes with clear commit messages.
3. Dry run before opening a PR:
   ```bash
   conda activate cva24v-env
   export NCBI_EMAIL="your@email"
   snakemake -s workflow/Snakefile -n -c 4
   ```
4. Verify `make test` works.
5. Open a pull request describing the change and any parameter updates.

## Code style
- Python scripts under `analysis/scripts` should use argparse and `--help` text.
- Prefer TSV for tables and FASTA for sequences.
- Avoid hard coded paths. Use flags and config.
