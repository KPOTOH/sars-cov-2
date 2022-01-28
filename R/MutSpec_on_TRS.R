setwd("~/COVID19/")
rm(list=ls(all=TRUE))

gis.trs_rs <- read.csv('data/gisaid_only_for_TRS-RS.csv')
gis.trs_rs$RefNuc = gis.trs_rs$parent_nucl
gis.trs_rs$AltNuc = gis.trs_rs$child_nucl

# read gisaid and filter not terminals
gis.full <- read.csv('data/gisaid_mutations_annotation.csv')
gis.full$RefNuc = gis.full$parent_nucl
gis.full$AltNuc = gis.full$child_nucl
gis.full <- subset(gis.full, RefNuc != '-' & AltNuc != '-')
gis.full$Pos <- gis.full$pos + 1
gis.full = subset(gis.full, ! startsWith(gis.full$child_node, '#'))

gis.utrs <- subset(gis.full, GenName %in% c('5UTR', '3UTR'))

annot <- read.csv('data/genes_annotation.csv')

trs_rs_pos <- c(29, 30, 31, 32, 33, 64, 72, 73, 74, 75, 25380, 25381, 25382, 
                25383, 25384, 25385, 25386, 25387, 25388, 25389, 26231, 26232, 
                26233, 26234, 26235, 26236, 26237, 26238, 27035, 27036, 27037, 
                27038, 27039, 27042, 27043, 27044, 27045, 27047, 27048, 29607, 
                29608, 29609, 29610, 29611, 29612, 29613, 29639, 29640, 29645, 29646)

ann.full = read.csv("data/ideal_table_of_sasha.csv")  # 1-based
ann.utrs <- subset(ann.full, GenName %in% c('5UTR', '3UTR'))

# ann.trans.black <- subset(ann.full, AaSub == 'S' & NucInCodon == 3 & Pos %in% trs_rs_pos)  # not true ff
# ann.untrans.black <- subset(ann.full, GenType == 'untranslated' & Pos %in% trs_rs_pos)

ann.black <- subset(ann.full, 
  (GenType == 'untranslated' & Pos %in% trs_rs_pos) | 
    (AaSub == 'S' & NucInCodon == 3 & Pos %in% trs_rs_pos)  # not 4f, but Syn and idx = 3
)

ann <- ann.black
df <- gis.trs_rs

mut_list = c('A>C','A>G','A>T','C>A','C>G','C>T','G>A','G>C','G>T','T>A','T>C','T>G')
mut_list = data.frame(mut_list)
names(mut_list)=c('NucSubst')
ann$MutExp <- 1
ann$NucSubst = paste(ann$RefNuc,ann$AltNuc,sep='>')
Exp12Comp = aggregate(ann$MutExp, by = list(ann$NucSubst), FUN = sum);
names(Exp12Comp) = c('NucSubst','ExpFr')

df$MutObs = 1
df$NucSubst = paste(df$RefNuc,df$AltNuc,sep='>')
Obs12Comp = aggregate(df$MutObs, by = list(df$NucSubst), FUN = sum);
names(Obs12Comp) = c('NucSubst','ObsFr')

MutSpec12Comp = merge(mut_list,Exp12Comp, by = 'NucSubst', all.x = TRUE)
MutSpec12Comp = merge(MutSpec12Comp,Obs12Comp, by = 'NucSubst',all.x = TRUE)
MutSpec12Comp[is.na(MutSpec12Comp)] <- 0
MutSpec12Comp$ObsToExp = MutSpec12Comp$ObsFr/MutSpec12Comp$ExpFr
MutSpec12Comp$ObsToExp[is.na(MutSpec12Comp$ObsToExp)] <- 0
MutSpec12Comp$MutSpec = MutSpec12Comp$ObsToExp/sum(MutSpec12Comp$ObsToExp)

barplot(MutSpec12Comp$MutSpec, names=MutSpec12Comp$NucSubst, cex.names = 0.7)
write.csv(x = MutSpec12Comp, 'data/MutSpecTables/07.MutSpec12_ForTRS_RS.csv')


# fisher.test(matrix(c(10.3, 4, 13.2, 3.5), nrow = 2))

