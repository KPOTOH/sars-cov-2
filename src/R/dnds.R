# insrall htslib from https://github.com/samtools/htslib (read errors...)
# BiocManager::install("Rhtslib")
# BiocManager::install("Rsamtools")
# install.packages("devtools", dependecies=TRUE)
# library(devtools); install_github("im3sanger/dndscv")

library("dndscv")

data("dataset_simbreast", package="dndscv")
dndsout = dndscv(mutations)

sel_cv = dndsout$sel_cv
print(head(sel_cv), digits = 3)

signif_genes = sel_cv[sel_cv$qglobal_cv==1, c("gene_name","qglobal_cv")]
rownames(signif_genes) = NULL
print(signif_genes)

path_cds_table = system.file("extdata", "BioMart_human_GRCh37_chr3_segment.txt", package = "dndscv", mustWork = TRUE)
path_genome_fasta = system.file("extdata", "chr3_segment.fa", package = "dndscv", mustWork = TRUE)

reftable = read.table(path_cds_table, header=1, sep="\t", stringsAsFactors=F)



setwd("~/COVID19")
dndscv::

mulal_path <- 'data/mulal_gisaid_2021-01-22.filtered.twice.fasta'
mulal = read.FASTA(mulal_path)
