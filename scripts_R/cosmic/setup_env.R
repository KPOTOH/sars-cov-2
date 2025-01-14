# on ubuntu
# sudo apt-get install -y libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev

renv::init(bioconductor = TRUE)

install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BiocStyle")
BiocManager::install("SomaticCancerAlterations")


devtools::install_github("nicolaroberts/hdp", build_vignettes = TRUE, dependencies = TRUE, force = TRUE)

# Suggests:
#   testthat,
#   SomaticCancerAlterations,
#   RColorBrewer,
#   knitr,
#   rmarkdown,
#   BiocStyle,
#   devtools


library("hdp")

hdp

