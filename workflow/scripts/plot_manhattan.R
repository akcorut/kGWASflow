## Logging
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## Load libraries
library(tidyverse)
library(qqman)

## Load input SAM alignment file of significant kmers (ignore SAM headers)
align_kmers_sam  <- read_delim(snakemake@input[[1]], comment = "@", 
                        col_names=FALSE, delim="\t")

## Select only columns with k-mer id, chromosome and chromosome position
kmers_align_tab <- align_kmers_sam %>%
    select(X1,X3,X4)

## Rename columns
colnames(kmers_align_tab) <- c("kmer", "chr", "bp")

## Split k-mer id column into k-mer number and p-value
kmers_align_tab <- kmers_align_tab %>%
    separate(kmer,  c("kmer", "p_value"), sep="_", remove=TRUE) %>%
    arrange(chr)

## Change column type
kmers_align_tab$p_value <- as.double(kmers_align_tab$p_value)   
kmers_align_tab$bp <- as.integer(kmers_align_tab$bp)
kmers_align_tab$chr <- as.factor(kmers_align_tab$chr)

## Get chromosome names in case they are in character format
## Adapted from: https://www.biostars.org/p/422781/
chrs <- unique(kmers_align_tab$chr)
kmers_align_tab$chr <- as.numeric(factor(kmers_align_tab$chr, levels = chrs))

## Get min p-value for the manhattan plot
min_pval <- min(kmers_align_tab$p_value)

print("Plotting...")
pdf(file = snakemake@output[["manhattan_plot"]], width = 12, height =6)
manhattan(kmers_align_tab, , chr = "chr", bp = "bp", p = "p_value", snp = "kmer",
            main = "Manhattan Plot", ylim=c(0, -log10(min_pval)+2), 
            cex = 1, cex.axis = 1, col = c("blue4"), 
            suggestiveline = F, genomewideline = F) 
            # chrlabs=chrs
dev.off()