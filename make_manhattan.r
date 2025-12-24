suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(dplyr))

# Check whether the provided arguments are valid
check_args_validity <- function(args) {
  # first check if the file exists
  if (!file.exists(args$input)) {
    print(paste0("ERROR: the input file, ", args$input, " does not exist"))
    stop()
  }
  # Next we need to ensure that if the user provided a pattern
  # that they also provided a directory and not a file
  is_file <- file.exists(args$input) && !dir.exists(args$input)
  if (is_file && args$pattern) {
    print("ERROR: A file and a pattern string were provided as inputs for the program. If you are providing a pattern as input, you need to provide a directory as the input.")
    stop()
  }
}

check_output_dir_exists <- function(output_dir) {
  if (!dir.exists(output_dir)) {
    print(paste0("The output directory, ", output_dir, " does not exist. Creating it now"))
    dir.create(args$output_dir)
  }
}

reformat <- function(dt, args) {
  # 1. Clean Chromosome column (remove 'chr', convert to numeric)
  dt[, CHR := as.numeric(gsub("chr", "", CHR, ignore.case = TRUE))]

  # 2. Calculate Cumulative Position (Vectorized - Instant)
  #    Get max BP per chromosome
  chr_max <- dt[, .(max_bp = max(POS, na.rm = TRUE)), by = CHR][order(CHR)]

  #    Calculate offset (lagged cumulative sum)
  chr_max[, offset := shift(cumsum(as.numeric(max_bp)), fill = 0, type = "lag")]

  #    Update main table by reference (fastest possible join)
  setkey(dt, CHR)
  setkey(chr_max, CHR)
  dt[chr_max, BPcum := POS + i.offset]

  #    Calculate axis labels
  axisdf <- dt[, .(center = (min(BPcum) + max(BPcum)) / 2), by = CHR]

  return(list(data = dt, axis = axisdf))
}

# If the user passes
gather_input_files <- function(input_path, pattern, pval_col) {
  if (dir.exists(input_path)) {
    input_files <- list.files(path = input_path, pattern = pattern, full.names = TRUE, recursive = T)
    print(paste0("Found ", length(input_files), " files from the analysis"))

    if (length(input_files) == 0) {
      print(paste0("There were not files found in the directory, ", input_path, " that had the pattern, ", pattern, ". Terminating program now."))
      stop()
    }
    print(paste0("Combining all ", length(input_files), " input files"))
    df <- input_files %>%
      map_dfr(~ fread(.)) %>%
      rename(c("POS" = "GENPOS", "PVAL" = "LOG10P", "CHR" = "CHROM"))
  } else if (!dir.exists(input_path) && file.exists(input_path)) {
    df <- fread(input_path) %>%
      rename(c("POS" = "GENPOS", "PVAL" = {{ pval_col }}, "CHR" = "CHROM"))
  }
  print(paste0("Read in ", nrow(df), " variants from the input"))
  return(df)
}


generate_manhattan <- function(snp_data, axisdf, args, colors) {

  snp_data[, is_annotate := ifelse(PVAL >= -log10(args$`annotate-threshold`), "yes", "no")]

  # Create label for annotation (Only create string for the few top hits)
  snp_data[is_annotate == "yes", SNP_ID := paste0(CHR, ":", POS)]

  manhattan <- ggplot(snp_data, aes(x = BPcum, y = PVAL)) +
    # Show all points
    geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
    scale_color_manual(values = rep(colors, 22)) +
    xlab("CHR") +
    # custom X axis:
    scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(-log10(snp_data$PVAL)) + 1)) +
    geom_hline(yintercept = -log10(args$bonferroni), color = "red", linetype = "dashed") +
    geom_hline(yintercept = -log10(args$suggestive), color = "blue", linetype = "dashed") + # remove space between plot area and x axis
    ggtitle(paste0("GWID:", args$`pheno-name`)) +
    # Custom the theme:
    theme_classic() +
    # Add label using ggrepel to avoid overlapping
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = args$`x-axis-text-size`),
      axis.text.y = element_text(size = args$`y-axis-text-size`),
      axis.title = element_text(size = args$`title-text-size`, face = "bold"),
      panel.grid.minor.x = element_blank()
    ) +
    geom_label_repel(data = snp_data[is_annotate == "yes"], aes(label = SNP_ID), size = 3)

  return(manhattan)
}

generate_qqplot <- function(snp_data, args) {

  raw_p <- 10^(-snp_data$PVAL)
  # Observed p-values
  observed <- sort(snp_data$PVAL)

  # Expected p-values (uniform distribution)
  expected <- sort(-log10(ppoints(length(observed))))

  qq_df <- data.table(observed = observed, expected = expected)

  # Lambda calculation
  chisq <- qchisq(1 - raw_p, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)

  p <- ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(size = 1) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    annotate("text", x = 1, y = max(observed) * 0.9, label = sprintf("λ = %.2f", lambda), size = 6) +
    labs(x = "Expected -log10(P)", y = "Observed -log10(P)") +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = args$`x-axis-text-size`),
      axis.text.y = element_text(size = args$`y-axis-text-size`),
    )

  return(p)
}

## Main script starts here: #########################
# Lets create a commandline parser
parser <- OptionParser()

parser <- add_option(parser, c("--input"), help = "Path to either the output file from the GWAS or to the directory that has the output files from the GWAS.", type = "character")

parser <- add_option(parser, c("--pheno-name"), help = "Name of the phenotype that the analysis was done for. This value will be used as the title in the output plots.", type = "character")

parser <- add_option(parser, c("--maf-threshold"), help = "Minor allele threshold so that SNPs below this value will be excluded from the plot.", type = "double", default = 0.01)

parser <- add_option(parser, c("--pattern"), help = "This is the pattern that will be used to glob files in the analysis. This argument should only be supplied if the user supplies a directory as the input.")

parser <- add_option(parser, c("--annotate-threshold"), help = "Minimum threshold to annotate SNPs on the Manhattan plot. These values should not be -log10 transformed. The default level is anything passing the Bonferroni threshold (5e-8).", type = "double", default = 5e-8)

parser <- add_option(parser, c("--bonferroni"), help = "Bonferroni threshold to use for the analysis. This value should not be -log10 transformed. (default= %default)", type = "double", default = 5e-8)

parser <- add_option(parser, c("--suggestive"), help = "Suggestive significance threshold. This value defaults to 1e-5 and should not be -log10 transformed.", type = "double", default = 1e-5)

parser <- add_option(parser, c("--output-dir"), help = "Path of a directory to write the output Manhattan and QQ plots to.", type = "character", default = "./test/")

parser <- add_option(parser, c("--pval-col"), help = "Name of the column that has p-values from the analysis.", type = "character", default = "PVal")

parser <- add_option(parser, c("--x-axis-text-size"), help = "Size of the x-axis title text.", type = "integer", default = 12)

parser <- add_option(parser, c("--y-axis-text-size"), help = "Size of the y-axis title text.", type = "integer", default = 12)

parser <- add_option(parser, c("--title-text-size"), help = "Size of the plot title text.", type = "integer", default = 14)

parser <- add_option(parser, c("--color-list"), help = "List of colors to be used in the Manhattan plots. Select only 2 colors and provide the colors as a comma-separated string with no spaces (Ex: 'red,blue'). Hexadecimal colors are also supported.", type = "character", default = "grey,skyblue")


args <- parse_args(parser)

check_args_validity(args)

check_output_dir_exists(args$`output-dir`)

## we can generate all of the config conditions
# split the color string
color_list <- unlist(strsplit(args$`color-list`, split = ","))

stopifnot(length(color_list) != 2, paste0("Expected only 2 colors to be provided to the program. Instead ", length(color_list), " colors were provided"))

df <- gather_input_files(args$input, args$pattern, args$`pval-col`)

df <- df %>% dplyr::filter(A1FREQ >= args$`maf-threshold`)

print(paste0(nrow(df), " variants passed the ", args$`maf-threshold`, " MAF threshold"))

# By this point we have changed the pvalue column to be PVAL
# reformatted_input = reformat(df, "PVAL")
reformatted_inputs <- reformat(df, args)

df_plot <- reformatted_inputs$data
axis_df <- reformatted_inputs$axis

print("Generating a manhattan plot from the GWAS results")
manhattan <- generate_manhattan(df_plot, axis_df, args, color_list)

ggsave(file.path(args$`output-dir`, paste0(args$`pheno-name`, "_gwas_manhattan_v2.png")), plot = manhattan, width = 12, height = 8)

print("Generating a QQ-plot from the GWAS results")
qqplot <- generate_qqplot(df_plot, args)

ggsave(file.path(args$`output-dir`, paste0(args$`pheno-name`, "_qqplot_v2.png")), plot = qqplot)
