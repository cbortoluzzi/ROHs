setwd("/path/to/input/file/from/step2/")

# Imports
library(qqman)

# Open all input files found in the directory
files <- list.files(pattern="*_input.txt")

# Plot genome-wide level of heterozygosity for each individual
for (input in files){
        file <- read.table(input,header=F,col.names =c("sample","chrom","start_bin","end_bin","nsites","nhet"))
        file$SNPcount <- round((10000/file$nsites)*file$nhet,2)
        pdf(paste0("heterozygosity_", file$sample,".pdf"), width=20)
        par(family="Times")
        plot <- manhattan(file, chr="chrom", bp="start_bin", p="SNPcount", snp="nhet", col="blue", suggestiveline=F, genomewideline=F, logp=F, type="h", ylab="Number of SNPs in bin", ylab="", ylim=c(0,200), cex.axis=1.0, cex.lab=1.0, main = unique(file$sample), cex.main=1.0)
        dev.off()
}

