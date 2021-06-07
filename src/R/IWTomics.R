## ----include=FALSE----------------------------------------------
library(knitr)
opts_chunk$set(fig.path='IWTomics-',concordance=TRUE,warning=TRUE,message=TRUE)
options(width=66)

# # ----installation,eval=FALSE,size="small"-----------------------
#  if (!requireNamespace("BiocManager", quietly=TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("IWTomics",dependencies=TRUE)

## ----loading,eval=TRUE,results='hide',warning=FALSE,message=FALSE,size="small"----
library(IWTomics)

# ----ETn_dataset,eval=TRUE,echo=TRUE,results='markup',size="small"----
data(ETn_example)
ETn_example

## ----ETn boxplot,eval=TRUE,echo=TRUE,fig.height=7,fig.width=10,label=ETn_boxplot,results='markup',fig.show='hide',size="small"----
plot(ETn_example,cex.main=2,cex.axis=1.2,cex.lab=1.2)

## ----ETn two samples test - mean,eval=TRUE,echo=TRUE,results='markup',size="small"----
ETn_test=IWTomicsTest(ETn_example,
                     id_region1='ETn_fixed',id_region2='Control')
adjusted_pval(ETn_test)

## ----ETn test plot,eval=TRUE,echo=TRUE,fig.height=10,fig.width=8,label=ETn_test_plot,results='markup',fig.show='hide',size="small"----
plotTest(ETn_test)
# plotSummary(ETn_test)

## ----ETn test plot scale 10kb,eval=TRUE,echo=TRUE,fig.height=10,fig.width=8,label=ETn_test_plot_10kb_scale,results='markup',fig.show='hide',size="small"----
adjusted_pval(ETn_test,scale_threshold=10)
plotTest(ETn_test,scale_threshold=10)

## ----datasets,eval=TRUE,echo=TRUE,results='markup',size="small"----
examples_path <- system.file("extdata",package="IWTomics")
datasets=read.table(file.path(examples_path,"datasets.txt"),
                    sep="\t",header=TRUE,stringsAsFactors=FALSE)
datasets

## ----feature datasets BED,eval=TRUE,echo=TRUE,results='markup',size="small"----
features_datasetsBED=
    read.table(file.path(examples_path,"features_datasetsBED.txt"),
               sep="\t",header=TRUE,stringsAsFactors=FALSE)
features_datasetsBED

## ----import features from BED,eval=TRUE,echo=TRUE,results='markup',size="small"----
startBED=proc.time()
regionsFeatures=IWTomicsData(datasets$regionFile,
    features_datasetsBED[,3:6],alignment='center',
    id_regions=datasets$id,name_regions=datasets$name,
    id_features=features_datasetsBED$id,
    name_features=features_datasetsBED$name,
    path=file.path(examples_path,'files'))
endBED=proc.time()
endBED-startBED

## ----feature datasets table,eval=TRUE,echo=TRUE,results='markup',size="small"----
features_datasetsTable=
    read.table(file.path(examples_path,"features_datasetsTable.txt"),
               sep="\t",header=TRUE,stringsAsFactors=FALSE)
features_datasetsTable

## ----import features from table,eval=TRUE,echo=TRUE,results='markup',size="small"----
startTable=proc.time()
regionsFeatures=IWTomicsData(datasets$regionFile,
    features_datasetsTable[,3:6],alignment='center',
    id_regions=datasets$id,name_regions=datasets$name,
    id_features=features_datasetsBED$id,
    name_features=features_datasetsBED$name,
    path=file.path(examples_path,'files'))
endTable=proc.time()
endTable-startTable

## ----IWTomicsData object,eval=TRUE,echo=TRUE,results='markup',size="small"----
regionsFeatures
regionsFeatures_subset=regionsFeatures[c('elem1','control'),'ftr1']
regionsFeatures_subset

## ----scatterplot,eval=TRUE,echo=TRUE,fig.height=5,fig.width=5,label=scatterplot,results='markup',fig.show='hide',size="small"----
plot(regionsFeatures,type='pairs')

## ----scatterplot smooth,eval=TRUE,echo=TRUE,fig.height=5,fig.width=5,label=scatterplot_smooth,results='markup',fig.show='hide',size="small"----
plot(regionsFeatures,type='pairsSmooth',col='violet')

## ----curves,eval=TRUE,echo=TRUE,fig.height=7,fig.width=10,label=curves,results='markup',fig.show='hide',size="small"----
plot(regionsFeatures,type='curves',
     N_regions=lengthRegions(regionsFeatures),
     id_features_subset='ftr1',cex.main=2,cex.axis=1.2,cex.lab=1.2)

## ----boxplot,eval=TRUE,echo=TRUE,fig.height=7,fig.width=10,label=boxplot,results='markup',fig.show='hide',size="small"----
plot(regionsFeatures,type='boxplot',
     id_features_subset='ftr2',cex.main=2,cex.axis=1.2,cex.lab=1.2)

## ----one sample test - mean,eval=TRUE,echo=TRUE,results='markup',size="small"----
result1=IWTomicsTest(regionsFeatures,mu=c(0,5),
                     id_region1='control',id_region2='',
                     id_features_subset=c('ftr1','ftr2'))
result1
adjusted_pval(result1)

## ----two samples test - mean,eval=TRUE,echo=TRUE,results='markup',size="small"----
result2_mean=IWTomicsTest(regionsFeatures,
                          id_region1=c('elem1','elem2','elem3'),
                          id_region2=c('control','control','control'))
result2_mean
adjusted_pval(result2_mean)

## ----two samples test - multiple quantiles,eval=TRUE,echo=TRUE,results='markup',size="small"----
result2_quantiles=IWTomicsTest(regionsFeatures,
                               id_region1=c('elem1','elem2','elem3'),
                               id_region2=c('control','control','control'),
                               id_features_subset='ftr1',
                               statistics='quantile',probs=c(0.25,0.75))
result2_quantiles
adjusted_pval(result2_quantiles)

## ----reproducibility,eval=TRUE,echo=TRUE,results='markup',size="small"----
set.seed(16)
result_rep1=IWTomicsTest(regionsFeatures,
                         id_region1='elem1',id_region2='control',
                         id_features_subset='ftr1')
adjusted_pval(result_rep1)

set.seed(16)
result_rep2=IWTomicsTest(regionsFeatures,
                         id_region1='elem1',id_region2='control',
                         id_features_subset='ftr1')
adjusted_pval(result_rep2)

identical(result_rep1,result_rep2)

## ----regionsFeatures_center,eval=TRUE,echo=TRUE,fig.height=7,fig.width=10,label=boxplot_diff_length,results='markup',fig.show='hide',size="small"----
data(regionsFeatures_center)
range(width(regions(regionsFeatures_center)))
plot_data=plot(regionsFeatures_center,type='boxplot',
               id_regions_subset=c('elem1','control'),
               id_features_subset='ftr1',size=TRUE)
plot_data$features_position_size

## ----not fully computable p-value curves,eval=TRUE,echo=TRUE,results='markup',warnings=TRUE,size="small"----
result_warning=IWTomicsTest(regionsFeatures_center,
                            id_region1='elem1',id_region2='control',
                            id_features_subset='ftr1')
adjusted_pval(result_warning)

## ----two samples test - mean - detailed plot,eval=TRUE,echo=TRUE,fig.height=10,fig.width=8,label=detailed_plot,results='markup',fig.show='hide',size="small"----
plotTest(result2_mean,alpha=0.05,id_features_subset='ftr1')

## ----two samples test - mean - summary plot feature,eval=TRUE,echo=TRUE,fig.height=2.5,fig.width=12,label=summary_plot_2_feature,results='markup',fig.show='hide',size="small"----
plotSummary(result2_mean,alpha=0.05,groupby='feature',align_lab='Center')

## ----one samples test - mean - summary plot feature,eval=TRUE,echo=TRUE,fig.height=2.5,fig.width=12,label=summary_plot_1_feature,results='markup',fig.show='hide',size="small"----
result1_mu1=IWTomicsTest(regionsFeatures,mu=1,
                         id_region1=c('elem1','elem2','elem3','control'))
plotSummary(result1_mu1,alpha=0.05,groupby='feature',align_lab='Center')

## ----one samples test - mean - summary plot location,eval=TRUE,echo=TRUE,fig.height=2,fig.width=12,label=summary_plot_1_location,results='markup',fig.show='hide',size="small"----
plotSummary(result1_mu1,alpha=0.05,groupby='test',align_lab='Center')

## ----many tests - mean - summary plot feature,eval=TRUE,echo=TRUE,fig.height=4,fig.width=12,label=summary_plot_many_feature,results='markup',fig.show='hide',size="small"----
result_many=IWTomicsTest(regionsFeatures,
                         id_region1=c('elem1','elem2','elem3',
                                      'elem1','elem1','elem2',
                                      'elem1','elem2','elem3','control'),
                         id_region2=c(rep('control',3),
                                      'elem2','elem3','elem3',
                                      rep('',4)),
                         id_features_subset='ftr1')
plotSummary(result_many,alpha=0.05,groupby='feature',
            align_lab='Center',gaps_tests=c(3,6),
            only_significant=TRUE)

## ----curve alignment scale,eval=TRUE,echo=TRUE,fig.height=5,fig.width=10,label=curves_scale,results='markup',fig.show='hide',size="small"----
data(regionsFeatures_scale)
lengthFeatures(regionsFeatures_scale[,'ftr1'])
features(regionsFeatures_scale)[['ftr1']][['elem1']][,1]
plot(regionsFeatures_scale,type='curves',
     N_regions=lengthRegions(regionsFeatures_scale),
     id_features_subset='ftr1')

## ----curve alignment scale - measurements on the same grid,eval=TRUE,echo=TRUE,fig.height=7,fig.width=10,label=curves_scale_same_grid,results='markup',fig.show='hide',size="small"----
regionsFeatures_scale_smooth=smooth(regionsFeatures_scale,type='locpoly')
lengthFeatures(regionsFeatures_scale_smooth[,'ftr1'])
features(regionsFeatures_scale_smooth)[['ftr1']][['elem1']][,1]
plot(regionsFeatures_scale_smooth,type='curves',
    N_regions=lengthRegions(regionsFeatures_scale_smooth),
     id_features_subset='ftr1')

## ----curve alignment scale - test,eval=TRUE,echo=TRUE,fig.height=2.5,fig.width=12,label=summary_plot_scale,results='markup',fig.show='hide',size="small"----
result=IWTomicsTest(regionsFeatures_scale_smooth,
                    id_region1=c('elem1','elem2','elem3'),
                    id_region2=c('control','control','control'),
                    id_features_subset='ftr1')
adjusted_pval(result)
plotSummary(result,groupby='feature')
# plotSummary(result)

## ----change resolution,eval=TRUE,echo=TRUE,results='markup',size="small"----
regionsFeatures
lengthFeatures(regionsFeatures)
regionsFeatures_smooth=smooth(regionsFeatures,type='locpoly',
                              resolution=c(4000,1000))
regionsFeatures_smooth
lengthFeatures(regionsFeatures_smooth)

## ----selection of relevant scale,eval=TRUE,echo=TRUE,fig.height=10,fig.width=8,label=relevant_scale,results='markup',fig.show='hide',size="small"----
result_scale=IWTomicsTest(regionsFeatures_center,
                          id_region1='elem1',id_region2='control',
                          id_features_subset='ftr1')
plotTest(result_scale,alpha=0.05,scale_threshold=8)

## ----setup,eval=TRUE,echo=TRUE,results='markup',size="small"----
sessionInfo()

