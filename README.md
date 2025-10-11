# pairwise-colocalization-smk

## Purpose

This pipeline performs pairwise colocalization analysis between genetic signals from different GWAS summary statistics using the [coloc](https://chr1swallace.github.io/coloc/) R package. It identifies pairs of signals that are close together, extracts the relevant regions from summary statistics, and runs colocalization for each pair.

## Overview

The workflow is managed by [Snakemake](https://snakemake.readthedocs.io/en/stable/) and is modular, allowing for scalable and reproducible analysis of many signal pairs.

## Inputs

- **Input signals file**: CSV file with columns:
  - `chromosome`: Chromosome number
  - `position`: Genomic position (base pair)
  - `gwas`: Name of the GWAS (should match summary stats file prefix)
  - (other columns allowed, but not required)
  - Example: `test/data/input_signals.csv`

- **Summary statistics files**: One file per GWAS, named `{gwas}.txt` and placed in the directory specified in `workflow/config.yaml` (default: `test/data/sumstats/`). Each file should be tab-delimited and contain at least:
  - `chromosome`
  - `position`
  - `rsid`
  - `beta`
  - `se`

## Outputs

- `outputs/pairs_to_test.csv`: List of all signal pairs to be tested
- `outputs/extracted_regions/`: Extracted region files for each pair and GWAS
- `outputs/coloc_results/`: Colocalization results for each pair
- `outputs/coloc_credible_sets/`: Credible set (if applicable) for each result
- `outputs/joined_credible_sets/`: For colocalized signals, credible sets are joined

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/en/stable/) (v6+ recommended)
- [R](https://www.r-project.org/) (v4+ recommended)
- R packages: `coloc`, `dplyr`, `readr`, `tidyr` 

## Usage

1. **Configure your input files**
   - Place your input signals CSV in `test/data/input_signals.csv` (or update the path in `workflow/config.yaml`)
   - Place your summary statistics files in `test/data/sumstats/` (or update the path in `workflow/config.yaml`)

2. **Edit configuration**
   - Edit `workflow/config.yaml` to set the correct paths and distance threshold (in base pairs)

3. **Run the pipeline**
   ```bash
   snakemake --cores 4 all
   ```
   (Adjust the number of cores as needed)

4. **View results**
   - Colocalization results for each pair will be in `test/outputs/coloc_results/`

## Pipeline Steps

1. **find_pairs**: Identifies all pairs of signals within the specified distance on the same chromosome, from different GWAS.
2. **extract_regions**: For each pair, extracts the relevant region from each GWAS summary stats file.
3. **run_coloc**: Runs the coloc R package on each pair of extracted regions and outputs the results.

## Customization

- You can adjust the distance threshold in `workflow/config.yaml`.
- You can modify the R scripts in `workflow/scripts/` to change how regions are extracted or how coloc is run.

## Troubleshooting

- Make sure all input files are present and paths are correct in the config.
- Check that your summary stats files have the required columns and are tab-delimited.
- If you see errors about missing R packages, install them in your R environment.

## Contact

For questions or issues, please contact the pipeline maintainer.
