
# I kept the original copy of Ryan's script intact (plot_res4.r)

suppressPackageStartupMessages(library(Haplin))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(CMplot))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(dplyr))

# Check whether the provided arguments are valid
check_args = function(args) {

}

# Lets create a commandline parser
parser = OptionParser()

parser = add_option(parser, c("--input"), help="Path to either the output file from the GWAS or to the directory that has the output files from the GWAS", type="character")

parser = add_option(parser, c("--phecode-name"), help="Name of the phenotype that the analysis was done for. This value will be used as the title in the output plots.", type="character")

parser = add_option(parser, c("--maf-threshold"), help="Minor allele threshold so that SNPs below this value will be excluded from the plot.", type="double", default=0.01)

parser = add_option(parser, c("--pattern"), help="This is the pattern that will be used to glob files in the analysis. This argument should only be supplied if the user supplies a directory as the input.")

parser = add_option(parser, c("--annotate-threshold"), help="Minumum threshold to annotate SNPs on the manhattan plot. These values should not be -log10 transformed. Default level is anything pass bonferroni (5e-8).", type="double", default=5e-8)

parser = add_option(parser, c("--bonferroni"), help="Bonferroni threshold to use for the analysis. This value should not be -log10 transformed. (default= %default)", type="double", default=5e-8)

parser = add_option(parser, c("--suggestive"), help="Suggestive significance threshold. This value defaults to 1e-5 and should not be -log10 transformed.", type="double", default=1e-5)

parser = add_option(parser, c("--output-dir"), help="Path of a directory to write the output Manhattan and QQ-plot to.", type="character", default="test.png")

# parser = add_option(parser)

parse_args(parser)

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

reformat = function(df, pval_col) {
    out = df %>% 
            select(c("ID", "CHR", "POS", {{ pval_col }})) %>%
            drop_na() %>%
            mutate("CHR" = as.numeric(gsub("chr", "", CHR))) 
    colnames(out) = c("SNP", "CHR", "BP", "P")
    out
}


MAF_THRESHOLD=0.001


# full_output_path = paste("./", , sep="/")
phecode_name = "Parkinsons"
input_dir = "/data100t1/home/james/agd250k/parkinsons/gwas/phers_results/step2_output"

# we need to merge the initial 22 chromosomes into 1 file 
input_files = list.files(path = input_dir, pattern = "*.regenie$", full.names=TRUE, recursive=T)
print(paste0("Found ", length(input_files), " files from the analysis"))

if (length(input_files) == 0) {
    quit()
}

combined_df = input_files %>%
                map_dfr(~fread(.)) %>% 
                rename(c("POS"="GENPOS", "PVAL"="LOG10P", "CHR"="CHROM")) 
                
combined_df = combined_df %>% dplyr::filter(A1FREQ >= MAF_THRESHOLD)

print(nrow(combined_df))

annolevel = -log10(5e-8)


df = combined_df %>% mutate("ID" = paste0(combined_df$CHR, ":", combined_df$POS))

reformatted_input = reformat(df, "PVAL")

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
  mutate( is_annotate=ifelse(P>annolevel, "yes", "no")) 

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
    geom_hline(yintercept=-log10(0.05/1000000), color="red", linetype="dashed") +  
    geom_hline(yintercept=-log10(1e-5), color="blue", linetype="dashed") +# remove space between plot area and x axis
    ggtitle(paste0("GWID:", phecode_name)) +
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

plots_dir = "/data100t1/home/james/agd250k/parkinsons/gwas/phers_results/plots"
if (!dir.exists(plots_dir)) {
    dir.create(plots_dir)
}


ggsave(file.path(plots_dir, "agd250k_parkinsons_phers_gwas_manhattan.png"), plot=manhattan, width=12, height=8)

reformatted_input = reformatted_input %>% mutate(OrigP=10**(-1*P)) %>% arrange(OrigP)

chisq <- qchisq(1 - reformatted_input$OrigP, 1)
lambda <- median(chisq) / qchisq(0.5, 1)

qqplot = ggplot(reformatted_input, aes(sample=P)) +
  stat_qq(distribution = stats::Uniform) + 
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

ggsave(file.path(plots_dir, "agd250k_parkinsons_qqplot_12_11_25.png"), plot=qqplot)


# plot_data <- reformatted_input %>%
#   arrange(OrigP) %>%   # Sort p-values from smallest to largest
#   mutate(
#     # Calculate Observed -log10 P-values
#     observed = -log10(OrigP),
    
#     # Calculate Expected -log10 P-values
#     # ppoints(n) generates n evenly spaced probabilities between 0 and 1
#     expected = -log10(ppoints(n()))
#   )

# lambda_value <- median(qchisq(1 - reformatted_input$OrigP, 1)) / qchisq(0.5, 1)
# lambda_label <- paste("Lambda =", round(lambda_value, 3))

# v2_qq = ggplot(plot_data, aes(x = expected, y = observed)) +
#   # Add the points
#   geom_point(alpha = 0.5, size = 2) +
  
#   # Add the null hypothesis line (y = x)
#   geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  
#   # Add labels and theme
#   labs(
#     title = "Parkinsons QQ Plot",
#     x = expression(Expected ~ -log[10](P)),
#     y = expression(Observed ~ -log[10](P)),
#     caption = lambda_label
#   ) +
#   theme_minimal() + 
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     panel.grid.minor = element_blank()
#   )

# ggsave(file.path(plots_dir, "agd250k_parkinsons_qqplot_v2_12_11_25.png"), plot=v2_qq)
# print(head(reformatted_input))
# # print(head(df))
# # This will pull out the value for the FWER or FDR
# # genomewide_avg = df[1, "PAdjCutoff"]
# # df$PVal is the list of p values we want to feed into the chisq
# # http://genometoolbox.blogspot.com/2014/08/how-to-calculate-genomic-inflation.html
# chisq <- qchisq(1 - reformatted_input$P, 1)
# avg <- median(chisq) / qchisq(0.5, 1)
# avg <- round(avg, digits = 3)
# pdf(output_qq_filename, 10, 10)
# pQQ(reformatted_input$P, main= paste(phecode_name,  "AFR QQ-plot"))
# text(0.5, 5, avg, cex = 1, font = 2)
# dev.off()

# gen.line = -log10(5e-8)
# # sug.line = NULL

# # df <- reformat(df)
# # print(head(df))
# minp = max(c(min(df$P), 0.05))
# pdf(output_manhattan_filename, 8, 8)



# manhattan(reformatted_input, main = paste(phecode_name,  "AFR Manhattan plot"), annotatePval = minp, genomewideline = gen.line, suggestiveline = -log10(1e-5), cex.lab=1.2, cex.axis=1.2)
# # manhattan(df, main = paste(phecode_name, "(", phecode, ")",  "Manhattan plot"), annotatePval = minp, genomewideline = -log10(genomewide_avg), ylim=c(0,8), suggestiveline = 0, cex.lab=1.2, cex.axis=1.2)

# # abline(h = -log10(1/100001), col = "blue", lty=2, lwd=3)  
# dev.off()   
    
# write phecode, phecode_name, and avg to a txt file 
# datastr = paste("", "\t", phecode_name, "\t", avg, "\n")
    # cat(datastr, file="genomic_inflation_factors.txt", append = TRUE)
    
# }
