
# I kept the original copy of Ryan's script intact (plot_res4.r)

suppressPackageStartupMessages(library(Haplin))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(CMplot))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(dplyr))

filter = function(df) {
    lastpos = 0
    lastchr = 0
    pos = 0
    out = NULL
    for(row in 1:dim(df)[1]) {
        chrom = df[row, 1]
        pos = pos + df[row, 3]
        if(chrom == lastchr) {
            if(pos - lastpos > 1) {
                lastchr = chrom
                lastpos = pos
                out = rbind(out, df[row,])
            }
        } else {
            lastchr = chrom
            lastpos = pos
            out = rbind(out, df[row,])
        }
    }
    out
}

# Check whether the provided arguments are valid
check_args_validity = function(args) {

  # first check if the file exists
  if (!file.exists(args$input)) {
    print(paste0("ERROR: the input file, ", args$input, " does not exist"))
    stop()
  }
  # Next we need to ensure that if the user provided a pattern 
  # that they also provided a directory and not a file
  is_file = file.exists(args$input) && !dir.exists(args$input)
  if (is_file &&  args$pattern) {
    print("ERROR: A file and a pattern string were provided as inputs for the program. If you are providing a pattern as input, you need to provide a directory as the input.")
    stop()
  } 
}

check_output_dir_exists = function(output_dir) {
  if (!dir.exists(output_dir)) {
    print(paste0("The output directory, ", output_dir, " does not exist. Creating it now"))
    dir.create(args$output_dir)
  }
}

reformat = function(df, pval_col) {
    df = df %>% mutate("ID" = paste0(df$CHR, ":", df$POS))
    out = df %>% 
            select(c("ID", "CHR", "POS", {{ pval_col }})) %>%
            drop_na() %>%
            mutate("CHR" = as.numeric(gsub("chr", "", CHR))) 
    colnames(out) = c("SNP", "CHR", "BP", "P")
    return(out)
}

# If the user passes 
gather_input_files = function(input_path, pattern, pval_col) {
  if dir.exists(input_path) {
    input_files = list.files(path = input_path, pattern = pattern, full.names=TRUE, recursive=T)
    print(paste0("Found ", length(input_files), " files from the analysis"))

    if (length(input_files) == 0) {
        print(paste0("There were not files found in the directory, ", input_path, " that had the pattern, ", pattern, ". Terminating program now."))
        stop()
    }
    print(paste0("Combining all ", length(input_file), " input files"))
    df = input_files %>%
                map_dfr(~fread(.)) %>% 
                rename(c("POS"="GENPOS", "PVAL"="LOG10P", "CHR"="CHROM")) 
  } else if (!dir.exists(input_path) && file.exists(input_path)) {
    df = fread(input_path) %>%
            rename(c("POS"="GENPOS", "PVAL"={{ pval_col }}, "CHR"="CHROM")) 
  }
  print(paste0("Read in" , nrow(combined_df), " variants from the input"))
  return(df)
}

generate_manhattan = function(args, df) {
  don <- reformatted_input %>% 
  
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(reformatted_input, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    mutate( is_annotate=ifelse(P>args$`annotate-threshold`, "yes", "no")) 

  axisdf = don %>%
    group_by(CHR) %>%
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

  manhattan = ggplot(don, aes(x=BPcum, y=P)) +
      # Show all points
      geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
      xlab("CHR") +
      # custom X axis:
      scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
      scale_y_continuous(limits=c(0,9), breaks=seq(0, 9, 1)) +   
      geom_hline(yintercept=-log10(args$bonferroni), color="red", linetype="dashed") +  
      geom_hline(yintercept=-log10(args$suggestive), color="blue", linetype="dashed") +# remove space between plot area and x axis
      ggtitle(paste0("GWID:", args$`phecode-name`)) +
      # Custom the theme:
      theme_classic() +
      # Add label using ggrepel to avoid overlapping
      geom_label_repel(data=subset(don, is_annotate=="yes"), aes(label=SNP), size=3) +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x= element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.minor.x = element_blank()
      )

  return(manhattan)
}

generate_qqplot = function(args, df) {
  reformatted_input = reformatted_input %>% mutate(OrigP=10**(-1*P)) %>% arrange(OrigP)

  chisq <- qchisq(1 - reformatted_input$OrigP, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)

  qqplot = ggplot(reformatted_input, aes(sample=P)) +
    stat_qq(distribution = stats::qunif) + 
    stat_qq_line() +
    xlab("observed") +
    ylab("expected") +
    annotate(
      geom = "text",
      x = -Inf,
      y = Inf,
      hjust = -0.15,
      vjust = 1 + 0.15 * 3,
      label = sprintf("λ = %.2f", lambda),
      size = 8
    ) +
    theme_classic() +
    theme(
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x= element_text(size=12),
      axis.text.y = element_text(size=12),
      axis.title=element_text(size=14,face="bold"),
      panel.grid.minor.x = element_blank()
    )
  return(qqplot)
}

# Lets create a commandline parser
parser = OptionParser()

parser = add_option(parser, c("--input"), help="Path to either the output file from the GWAS or to the directory that has the output files from the GWAS", type="character")

parser = add_option(parser, c("--pheno-name"), help="Name of the phenotype that the analysis was done for. This value will be used as the title in the output plots.", type="character")

parser = add_option(parser, c("--maf-threshold"), help="Minor allele threshold so that SNPs below this value will be excluded from the plot.", type="double", default=0.01)

parser = add_option(parser, c("--pattern"), help="This is the pattern that will be used to glob files in the analysis. This argument should only be supplied if the user supplies a directory as the input.")

parser = add_option(parser, c("--annotate-threshold"), help="Minumum threshold to annotate SNPs on the manhattan plot. These values should not be -log10 transformed. Default level is anything pass bonferroni (5e-8).", type="double", default=5e-8)

parser = add_option(parser, c("--bonferroni"), help="Bonferroni threshold to use for the analysis. This value should not be -log10 transformed. (default= %default)", type="double", default=5e-8)

parser = add_option(parser, c("--suggestive"), help="Suggestive significance threshold. This value defaults to 1e-5 and should not be -log10 transformed.", type="double", default=1e-5)

parser = add_option(parser, c("--output-dir"), help="Path of a directory to write the output Manhattan and QQ-plot to.", type="character", default="test.png")

parser = add_option(parser, c("--pval-col"), help="Name of the column that has pvalues from the analysis", type="character", default="PVal")

# parser = add_option(parser)

args = parse_args(parser)

check_args_validity(args)

check_output_dir_exists(args$`output-dir`)

df = gather_input_files(args$input, args$pattern, args$`pval-col`)

combined_df = df %>% dplyr::filter(A1FREQ >= args$`maf-threshold`)

print(paste0(nrow(combined_df), " variants that passed the ", args$`maf-threshold`, " MAF threshold"))

# By this point we have changed the pvalue column to be PVAL
reformatted_input = reformat(df, "PVAL")

manhattan = generate_manhattan(args, reformatted_input)

ggsave(file.path(args$`output-dir`, paste0(args$`pheno-name`, "_gwas_manhattan.png")), plot=manhattan, width=12, height=8)

qqplot = generate_qqplot(args, reformatted_input)

ggsave(file.path(args$`output-dir`, paste0(args$`pheno-name`,"agd250k_parkinsons_qqplot_12_11_25.png")), plot=qqplot)







