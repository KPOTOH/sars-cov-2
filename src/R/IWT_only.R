# ## ----include=FALSE----------------------------------------------
# library(knitr)
# opts_chunk$set(fig.path='IWTomics-',concordance=TRUE,warning=TRUE,message=TRUE)
# options(width=66)

library(IWTomics)
rm(list=ls(all=TRUE))
setwd("~/COVID19/data/_iwt_dataset20_10_10000_30000_lvl-2_mode-more/")

datasets=read.table("datasets.txt", sep="\t",header=TRUE,stringsAsFactors=FALSE)
features_datasetsTable=read.table("features_datasetsTable.txt", sep="\t", header=TRUE,stringsAsFactors=FALSE)
regionsFeatures=IWTomicsData(datasets$regionFile,
                             features_datasetsTable[,3:4],
                             # alignment='center',
                             id_regions=datasets$id,
                             name_regions=datasets$name,
                             id_features=features_datasetsTable$id,
                             name_features=features_datasetsTable$name,
                             # path=file.path(examples_path,'files'),
                             path='files',
)

plot_data = plot(regionsFeatures,type='boxplot',
                 id_features_subset='ftr1',
                 cex.main=2,cex.axis=1.2,cex.lab=1.2, ask=FALSE)

plot_data$features_plot  # this is quantiles


result2_quantiles=IWTomicsTest(regionsFeatures,
                               id_region1=c('elem1'),
                               id_region2=c('control'),
                               id_features_subset='ftr1',
                               statistics='quantile',probs=c(0.25, 0.5, 0.75))
adjusted_pval(result2_quantiles)




result2_mean=IWTomicsTest(regionsFeatures,
                          id_region1=c('elem1'),
                          id_region2=c('control'),
                          id_features_subset='ftr1')
adjusted_pval(result2_mean)

# plotTest(result2_mean)
