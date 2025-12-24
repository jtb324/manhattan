## make_mahattan script

This repository contains the script that I use to generate manhattan plots and to generate qqplots. This script is larger designed to work with REGENIE at the moment but *hopefully* will support more programs as time goes on.

## Installation

To run this script, you need to have R installed along with the following R packages:

- `optparse`
- `data.table`
- `ggplot2`
- `tidyverse`
- `ggrepel`
- `dplyr`

You can install these dependencies by running the following command in your R console:

```r
install.packages(c("optparse", "data.table", "tidyverse", "ggrepel"))
```

*Note: `dplyr` and `ggplot2` are included in `tidyverse`.*

Alternatively, you can use the provided `environment.yml` file to create a Conda environment named "manhatten_env" with all dependencies installed:

```bash
# Using conda
conda env create -f environment.yml
conda activate manhattan_env

# Or using mamba (recommended for speed)
mamba env create -f environment.yml
mamba activate manhattan_env
```

## Usage

The script is run from the command line using `Rscript`.

### Input Data Requirements

The input data (single file) should contain the following columns:

- `CHROM`: Chromosome number (e.g., "1", "chr1")
- `GENPOS`: Genomic position
- `A1FREQ`: Allele frequency
- A p-value column (name specified by `--pval-col`, default is "PVal")

### Example Command

Here is an example of how to run the script for a single GWAS result file:

```bash
Rscript make_manhattan.r \
  --input "path/to/gwas_results.txt" \
  --pheno-name "TraitName" \
  --output-dir "output_plots" \
  --pval-col "P_VALUE" \
  --maf-threshold 0.01 \
  --color-list "grey,skyblue"
```

### Arguments

- `--input`: Path to the GWAS results file or directory.
- `--pheno-name`: Name of the phenotype (used for plot titles and filenames).
- `--output-dir`: Directory where the output plots will be saved.
- `--pval-col`: Name of the column containing p-values (default: "PVal").
- `--maf-threshold`: Minor Allele Frequency threshold for filtering (default: 0.01).
- `--bonferroni`: Bonferroni significance threshold (default: 5e-8).
- `--suggestive`: Suggestive significance threshold (default: 1e-5).
- `--pattern`: (Optional) File pattern to match if input is a directory.
- `--x-axis-text-size`: Size of the x-axis title text (default: 12).
- `--y-axis-text-size`: Size of the y-axis title text (default: 12).
- `--title-text-size`: Size of the plot title text (default: 14).
- `--color-list`: List of colors for the Manhattan plot, provided as a comma-separated string (default: "grey,skyblue").
