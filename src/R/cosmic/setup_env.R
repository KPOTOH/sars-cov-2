# on ubuntu
# sudo apt-get install libcurl4-openssl-dev
# sudo apt-get install libxml2-dev
# sudo apt-get install libgit2-dev


install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BiocStyle")
BiocManager::install("SomaticCancerAlterations")


devtools::install_github("nicolaroberts/hdp", build_vignettes = TRUE, dependencies = TRUE)

# Suggests:
#   testthat,
#   SomaticCancerAlterations,
#   RColorBrewer,
#   knitr,
#   rmarkdown,
#   BiocStyle,
#   devtools