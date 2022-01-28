library(ggplot2)
library(lubridate)
library(dplyr)
library(plotly)

ideal_table = read.csv("D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data_obtained/ideal_table.csv")
#parent nuc == ref nuc
#gis = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/for_mut_spec_gis.csv')
gis = read.csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/gisaid_mutations_annotation.csv')
gis$pos = gis$pos+1
ideal_table[[1]] <- NULL
ideal_table[[1]] <- NULL
ideal_table$IsStem = ifelse(is.na(ideal_table$IsStem), yes = 0, no = 1)
ideal_table$Quadr = ifelse(is.na(ideal_table$Quadr), yes = 0, no = 1)

Stem_pos = ideal_table[ideal_table$IsStem == 1,]$Pos
Stem_pos = unique(Stem_pos)
length(Stem_pos)
Qudr_pos = ideal_table[ideal_table$Quadr == 1,]$Pos
Qudr_pos = unique(Qudr_pos)

gis_aa_change = gis[gis$parent_nucl!='-' & gis$child_nucl!='-',]
gis_aa_change = gis_aa_change[!is.na(gis_aa_change$pos),]
gis_aa_change[gis_aa_change$child_codon == 'CTT' | gis_aa_change$child_codon == 'CTA' | gis_aa_change$child_codon == 'CTG' | gis_aa_change$child_codon == 'CTC',]$child_aa = 'L_CT'
gis_aa_change[gis_aa_change$child_codon == 'TTA' | gis_aa_change$child_codon == 'TTG',]$child_aa = 'L_TT'

gis_aa_change[gis_aa_change$child_codon == 'TCT' | gis_aa_change$child_codon == 'TCA' | gis_aa_change$child_codon == 'TCG' | gis_aa_change$child_codon == 'TCC',]$child_aa = 'S_TC'
gis_aa_change[gis_aa_change$child_codon == 'AGT' | gis_aa_change$child_codon == 'AGC',]$child_aa = 'S_AG'

gis_aa_change[gis_aa_change$child_codon == 'CGC' | gis_aa_change$child_codon == 'CGA' | gis_aa_change$child_codon == 'CGT' | gis_aa_change$child_codon == 'CGG',]$child_aa = 'R_CG'
gis_aa_change[gis_aa_change$child_codon == 'AGG' | gis_aa_change$child_codon == 'AGA',]$child_aa = 'R_AG'

gis_aa_change[gis_aa_change$child_codon == 'TAG' | gis_aa_change$child_codon == 'TAA',]$child_aa = '*_TA'
gis_aa_change[gis_aa_change$child_codon == 'TGA',]$child_aa = '*_TG'

gis_aa_change$AaSub = ifelse(gis_aa_change$child_aa == gis_aa_change$parent_aa, 'S', 'NS')
gis_aa_change[gis_aa_change$GenType == 'untranslated',]$AaSub = 'NS'
gis_aa_change$MutObs = 1
gis_aa_change$NucSubst = paste(gis_aa_change$parent_nucl, gis_aa_change$child_nucl, sep = '>')
gis_syn_stem = gis_aa_change[gis_aa_change$pos %in% Stem_pos & gis_aa_change$AaSub == 'S',]
gis_all_stem = gis_aa_change[gis_aa_change$pos %in% Stem_pos,]
gis_syn_qudr = gis_aa_change[gis_aa_change$pos %in% Qudr_pos & gis_aa_change$AaSub == 'S',]
gis_all_qudr = gis_aa_change[gis_aa_change$pos %in% Qudr_pos,]

mut_list = c('A>C','A>G','A>T','C>A','C>G','C>T','G>A','G>C','G>T','T>A','T>C','T>G')
mut_list = data.frame(mut_list)
names(mut_list)=c('NucSubst')

#stem syn
ann = ideal_table[ideal_table$IsStem == 1 & ideal_table$AaSub == 'S',]
ann$MutExp = 1
ann$NucSubst = paste(ann$RefNuc, ann$AltNuc, sep = '>')
Exp12Comp = aggregate(ann$MutExp, by = list(ann$NucSubst), FUN = sum);
names(Exp12Comp) = c('NucSubst','ExpFr')

Obs12Comp = aggregate(gis_syn_stem$MutObs, by = list(gis_syn_stem$NucSubst), FUN = sum);
names(Obs12Comp) = c('NucSubst','ObsFr')

MutSpec12Comp = merge(mut_list,Exp12Comp, by = 'NucSubst', all.x = TRUE)
MutSpec12Comp = merge(MutSpec12Comp,Obs12Comp, by = 'NucSubst',all.x = TRUE)
MutSpec12Comp[is.na(MutSpec12Comp)] <- 0
MutSpec12Comp$ObsToExp = MutSpec12Comp$ObsFr/MutSpec12Comp$ExpFr
MutSpec12Comp$ObsToExp[is.na(MutSpec12Comp$ObsToExp)] <- 0
MutSpec12Comp$MutSpec = MutSpec12Comp$ObsToExp/sum(MutSpec12Comp$ObsToExp)

write.csv(x = MutSpec12Comp, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Stem/data/18.MutSpec_12_stem_syn.csv')

png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Stem/figures/18.MutSpec_12_stem_syn.png')
barplot(MutSpec12Comp$MutSpec, names =  MutSpec12Comp$NucSubst, main = 'Syn MutSpec 12 components for stems', cex.names = 0.7)
dev.off()

#stem
ann = ideal_table[ideal_table$IsStem == 1,]
ann$MutExp = 1
ann$NucSubst = paste(ann$RefNuc, ann$AltNuc, sep = '>')
Exp12Comp = aggregate(ann$MutExp, by = list(ann$NucSubst), FUN = sum);
names(Exp12Comp) = c('NucSubst','ExpFr')

Obs12Comp = aggregate(gis_all_stem$MutObs, by = list(gis_all_stem$NucSubst), FUN = sum);
names(Obs12Comp) = c('NucSubst','ObsFr')

MutSpec12Comp = merge(mut_list,Exp12Comp, by = 'NucSubst', all.x = TRUE)
MutSpec12Comp = merge(MutSpec12Comp,Obs12Comp, by = 'NucSubst',all.x = TRUE)
MutSpec12Comp[is.na(MutSpec12Comp)] <- 0
MutSpec12Comp$ObsToExp = MutSpec12Comp$ObsFr/MutSpec12Comp$ExpFr
MutSpec12Comp$ObsToExp[is.na(MutSpec12Comp$ObsToExp)] <- 0
MutSpec12Comp$MutSpec = MutSpec12Comp$ObsToExp/sum(MutSpec12Comp$ObsToExp)

write.csv(x = MutSpec12Comp, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Stem/data/18.MutSpec_12_stem_all.csv')

png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Stem/figures/18.MutSpec_12_stem_all.png')
barplot(MutSpec12Comp$MutSpec, names =  MutSpec12Comp$NucSubst, main = 'MutSpec 12 components for stems', cex.names = 0.7)
dev.off()

#quadr syn

ann = ideal_table[ideal_table$Quadr == 1 & ideal_table$AaSub == 'S',]
ann$MutExp = 1
ann$NucSubst = paste(ann$RefNuc, ann$AltNuc, sep = '>')
Exp12Comp = aggregate(ann$MutExp, by = list(ann$NucSubst), FUN = sum);
names(Exp12Comp) = c('NucSubst','ExpFr')

Obs12Comp = aggregate(gis_syn_qudr$MutObs, by = list(gis_syn_qudr$NucSubst), FUN = sum);
names(Obs12Comp) = c('NucSubst','ObsFr')

MutSpec12Comp = merge(mut_list,Exp12Comp, by = 'NucSubst', all.x = TRUE)
MutSpec12Comp = merge(MutSpec12Comp,Obs12Comp, by = 'NucSubst',all.x = TRUE)
MutSpec12Comp[is.na(MutSpec12Comp)] <- 0
MutSpec12Comp$ObsToExp = MutSpec12Comp$ObsFr/MutSpec12Comp$ExpFr
MutSpec12Comp$ObsToExp[is.na(MutSpec12Comp$ObsToExp)] <- 0
MutSpec12Comp$MutSpec = MutSpec12Comp$ObsToExp/sum(MutSpec12Comp$ObsToExp)

write.csv(x = MutSpec12Comp, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Stem/data/18.MutSpec_12_quadr_syn.csv')

png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Stem/figures/18.MutSpec_12_quadr_syn.png')
barplot(MutSpec12Comp$MutSpec, names =  MutSpec12Comp$NucSubst, main = 'Syn MutSpec 12 components for quadruplexes', cex.names = 0.7)
dev.off()

#quadr
ann = ideal_table[ideal_table$Quadr == 1,]
ann$MutExp = 1
ann$NucSubst = paste(ann$RefNuc, ann$AltNuc, sep = '>')
Exp12Comp = aggregate(ann$MutExp, by = list(ann$NucSubst), FUN = sum);
names(Exp12Comp) = c('NucSubst','ExpFr')

Obs12Comp = aggregate(gis_all_qudr$MutObs, by = list(gis_all_qudr$NucSubst), FUN = sum);
names(Obs12Comp) = c('NucSubst','ObsFr')

MutSpec12Comp = merge(mut_list,Exp12Comp, by = 'NucSubst', all.x = TRUE)
MutSpec12Comp = merge(MutSpec12Comp,Obs12Comp, by = 'NucSubst',all.x = TRUE)
MutSpec12Comp[is.na(MutSpec12Comp)] <- 0
MutSpec12Comp$ObsToExp = MutSpec12Comp$ObsFr/MutSpec12Comp$ExpFr
MutSpec12Comp$ObsToExp[is.na(MutSpec12Comp$ObsToExp)] <- 0
MutSpec12Comp$MutSpec = MutSpec12Comp$ObsToExp/sum(MutSpec12Comp$ObsToExp)

write.csv(x = MutSpec12Comp, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Stem/data/18.MutSpec_12_quadr_all.csv')

png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Stem/figures/18.MutSpec_12_quadr_all.png')
barplot(MutSpec12Comp$MutSpec, names =  MutSpec12Comp$NucSubst, main = 'MutSpec 12 components for quadruplexes', cex.names = 0.7)
dev.off()

#not stem not quadr
`%notin%` <- Negate(`%in%`)
gis_syn_notstem = gis_aa_change[gis_aa_change$pos %notin% Stem_pos & gis_aa_change$AaSub == 'S',]
gis_all_notstem = gis_aa_change[gis_aa_change$pos %notin% Stem_pos,]
gis_syn_notqudr = gis_aa_change[gis_aa_change$pos %notin% Qudr_pos & gis_aa_change$AaSub == 'S',]
gis_all_notqudr = gis_aa_change[gis_aa_change$pos %notin% Qudr_pos,]

#stem syn
ann = ideal_table[ideal_table$IsStem == 0 & ideal_table$AaSub == 'S',]
ann$MutExp = 1
ann$NucSubst = paste(ann$RefNuc, ann$AltNuc, sep = '>')
Exp12Comp = aggregate(ann$MutExp, by = list(ann$NucSubst), FUN = sum);
names(Exp12Comp) = c('NucSubst','ExpFr')

Obs12Comp = aggregate(gis_syn_notstem$MutObs, by = list(gis_syn_notstem$NucSubst), FUN = sum);
names(Obs12Comp) = c('NucSubst','ObsFr')

MutSpec12Comp = merge(mut_list,Exp12Comp, by = 'NucSubst', all.x = TRUE)
MutSpec12Comp = merge(MutSpec12Comp,Obs12Comp, by = 'NucSubst',all.x = TRUE)
MutSpec12Comp[is.na(MutSpec12Comp)] <- 0
MutSpec12Comp$ObsToExp = MutSpec12Comp$ObsFr/MutSpec12Comp$ExpFr
MutSpec12Comp$ObsToExp[is.na(MutSpec12Comp$ObsToExp)] <- 0
MutSpec12Comp$MutSpec = MutSpec12Comp$ObsToExp/sum(MutSpec12Comp$ObsToExp)

write.csv(x = MutSpec12Comp, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Stem/data/18.MutSpec_12_notstem_syn.csv')

png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Stem/figures/18.MutSpec_12_notstem_syn.png')
barplot(MutSpec12Comp$MutSpec, names =  MutSpec12Comp$NucSubst, main = 'Syn MutSpec 12 components for not stems', cex.names = 0.7)
dev.off()

#stem
ann = ideal_table[ideal_table$IsStem == 0,]
ann$MutExp = 1
ann$NucSubst = paste(ann$RefNuc, ann$AltNuc, sep = '>')
Exp12Comp = aggregate(ann$MutExp, by = list(ann$NucSubst), FUN = sum);
names(Exp12Comp) = c('NucSubst','ExpFr')

Obs12Comp = aggregate(gis_all_notstem$MutObs, by = list(gis_all_notstem$NucSubst), FUN = sum);
names(Obs12Comp) = c('NucSubst','ObsFr')

MutSpec12Comp = merge(mut_list,Exp12Comp, by = 'NucSubst', all.x = TRUE)
MutSpec12Comp = merge(MutSpec12Comp,Obs12Comp, by = 'NucSubst',all.x = TRUE)
MutSpec12Comp[is.na(MutSpec12Comp)] <- 0
MutSpec12Comp$ObsToExp = MutSpec12Comp$ObsFr/MutSpec12Comp$ExpFr
MutSpec12Comp$ObsToExp[is.na(MutSpec12Comp$ObsToExp)] <- 0
MutSpec12Comp$MutSpec = MutSpec12Comp$ObsToExp/sum(MutSpec12Comp$ObsToExp)

write.csv(x = MutSpec12Comp, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Stem/data/18.MutSpec_12_notstem_all.csv')

png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Stem/figures/18.MutSpec_12_notstem_all.png')
barplot(MutSpec12Comp$MutSpec, names =  MutSpec12Comp$NucSubst, main = 'MutSpec 12 components for no stems', cex.names = 0.7)
dev.off()

#quadr syn

ann = ideal_table[ideal_table$Quadr == 0 & ideal_table$AaSub == 'S',]
ann$MutExp = 1
ann$NucSubst = paste(ann$RefNuc, ann$AltNuc, sep = '>')
Exp12Comp = aggregate(ann$MutExp, by = list(ann$NucSubst), FUN = sum);
names(Exp12Comp) = c('NucSubst','ExpFr')

Obs12Comp = aggregate(gis_syn_notqudr$MutObs, by = list(gis_syn_notqudr$NucSubst), FUN = sum);
names(Obs12Comp) = c('NucSubst','ObsFr')

MutSpec12Comp = merge(mut_list,Exp12Comp, by = 'NucSubst', all.x = TRUE)
MutSpec12Comp = merge(MutSpec12Comp,Obs12Comp, by = 'NucSubst',all.x = TRUE)
MutSpec12Comp[is.na(MutSpec12Comp)] <- 0
MutSpec12Comp$ObsToExp = MutSpec12Comp$ObsFr/MutSpec12Comp$ExpFr
MutSpec12Comp$ObsToExp[is.na(MutSpec12Comp$ObsToExp)] <- 0
MutSpec12Comp$MutSpec = MutSpec12Comp$ObsToExp/sum(MutSpec12Comp$ObsToExp)

write.csv(x = MutSpec12Comp, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Stem/data/18.MutSpec_12_notquadr_syn.csv')

png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Stem/figures/18.MutSpec_12_notquadr_syn.png')
barplot(MutSpec12Comp$MutSpec, names =  MutSpec12Comp$NucSubst, main = 'Syn MutSpec 12 components for no quadruplexes', cex.names = 0.7)
dev.off()

#quadr
ann = ideal_table[ideal_table$Quadr == 0,]
ann$MutExp = 1
ann$NucSubst = paste(ann$RefNuc, ann$AltNuc, sep = '>')
Exp12Comp = aggregate(ann$MutExp, by = list(ann$NucSubst), FUN = sum);
names(Exp12Comp) = c('NucSubst','ExpFr')

Obs12Comp = aggregate(gis_all_notqudr$MutObs, by = list(gis_all_notqudr$NucSubst), FUN = sum);
names(Obs12Comp) = c('NucSubst','ObsFr')

MutSpec12Comp = merge(mut_list,Exp12Comp, by = 'NucSubst', all.x = TRUE)
MutSpec12Comp = merge(MutSpec12Comp,Obs12Comp, by = 'NucSubst',all.x = TRUE)
MutSpec12Comp[is.na(MutSpec12Comp)] <- 0
MutSpec12Comp$ObsToExp = MutSpec12Comp$ObsFr/MutSpec12Comp$ExpFr
MutSpec12Comp$ObsToExp[is.na(MutSpec12Comp$ObsToExp)] <- 0
MutSpec12Comp$MutSpec = MutSpec12Comp$ObsToExp/sum(MutSpec12Comp$ObsToExp)

write.csv(x = MutSpec12Comp, 'D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Stem/data/18.MutSpec_12_notquadr_all.csv')

png('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/R_figure/Stem/figures/18.MutSpec_12_notquadr_all.png')
barplot(MutSpec12Comp$MutSpec, names =  MutSpec12Comp$NucSubst, main = 'MutSpec 12 components for no quadruplexes', cex.names = 0.7)
dev.off()