library(ggplot2)
library(lubridate)
library(dplyr)
library(plotly)

setwd("~/COVID19/")
rm(list=ls(all=TRUE))

draw_mutspec = function(df, ann, n_mutspec = 12, mut_type = 'S', spec_type = 'general', 
                        region_type = 'translated', label = 'AllPos'){
  ann = as.data.frame(ann)
  df = as.data.frame(df)
  if (region_type == 'translated'){
    df = df[df$GenType == 'translated',]
    ann = ann[ann$GenType == 'translated',]
  }else if(region_type == 'untranslated'){
    df = df[df$GenType == 'untranslated',]
    ann = ann[ann$GenType == 'untranslated',]
  }
  if (mut_type == 'S'){
    df = df[df$AaSub == 'S',]
    ann = ann[ann$AaSub == 'S',]
  }else if(mut_type == 'NS'){
    df = df[df$AaSub == 'NS',]
    ann = ann[ann$AaSub == 'NS',]
  }
  if (spec_type == 'general'){
    if(n_mutspec == 12){
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
      
      title <- sprintf('07.MutSpec12_For%s', label)
      PATH_TO_MUTSPEC_TABLE <- sprintf('data/MutSpecTables/%s.csv', title)
      write.csv(x = MutSpec12Comp, PATH_TO_MUTSPEC_TABLE)
      
      png(sprintf('visual_data/MutSpecFigures/%s.png', title))
      barplot(MutSpec12Comp$MutSpec, names=MutSpec12Comp$NucSubst, 
              main=sprintf('MutSpec 12 components for %s', label), cex.names = 0.7)
      dev.off()
    }
  }else if(spec_type == 'genes'){
    if(n_mutspec == 12){
      for (gen in unique(df$GenName)){
        gen_df = df[df$GenName == gen,]
        gen_ann = ann[ann$GenName == gen,]
        mut_list = c('A>C','A>G','A>T','C>A','C>G','C>T','G>A','G>C','G>T','T>A','T>C','T>G')
        mut_list = data.frame(mut_list)
        names(mut_list)=c('NucSubst')
        gen_ann$MutExp = 1
        gen_ann$NucSubst = paste(gen_ann$RefNuc,gen_ann$AltNuc,sep='>')
        Exp12Comp = aggregate(gen_ann$MutExp, by = list(gen_ann$NucSubst), FUN = sum);
        names(Exp12Comp) = c('NucSubst','ExpFr')
        
        gen_df$MutObs = 1
        gen_df$NucSubst = paste(gen_df$RefNuc,gen_df$AltNuc,sep='>')
        Obs12Comp = aggregate(gen_df$MutObs, by = list(gen_df$NucSubst), FUN = sum);
        names(Obs12Comp) = c('NucSubst','ObsFr')
        
        MutSpec12Comp = merge(mut_list,Exp12Comp, by = 'NucSubst', all.x = TRUE)
        MutSpec12Comp = merge(MutSpec12Comp,Obs12Comp, by = 'NucSubst',all.x = TRUE)
        MutSpec12Comp[is.na(MutSpec12Comp)] <- 0
        MutSpec12Comp$ObsToExp = MutSpec12Comp$ObsFr/MutSpec12Comp$ExpFr
        MutSpec12Comp$ObsToExp[is.na(MutSpec12Comp$ObsToExp)] <- 0
        MutSpec12Comp$MutSpec = MutSpec12Comp$ObsToExp/sum(MutSpec12Comp$ObsToExp)
        write.csv(x = MutSpec12Comp, sprintf('data/MutSpecTables/07.MutSpec12_%s_%s.csv', label, gen))
        
        png(sprintf('visual_data/MutSpecFigures/07.MutSpec12_%s_%s.png', label, gen))
        barplot(MutSpec12Comp$MutSpec, names=MutSpec12Comp$NucSubst, 
                main=sprintf('MutSpec 12 components per Gen %s in %s', gen, label), cex.names = 0.7)
        dev.off()
      }
    }
  }
}

# read data
ann = read.csv("data/ideal_table_of_sasha.csv")  # 1-based
ffcodons = c('TCA', 'TCT', 'TCG', 'TCC',
             'CTA', 'CTT', 'CTG', 'CTC', 
             'CCA', 'CCT', 'CCG', 'CCC', 
             'CGA', 'CGT', 'CGG', 'CGC', 
             'ACA', 'ACT', 'ACG', 'ACC', 
             'GTA', 'GTT', 'GTG', 'GTC', 
             'GCA', 'GCT', 'GCG', 'GCC', 
             'GGA', 'GGT', 'GGG', 'GGC')

ann.ff = subset(ann, NucInCodon == 3 & AaSub == 'S' & RefCodon %in% ffcodons)

# gis = read.csv('data/gisaid_mutations_annotation.csv')  # 0-based
gis = read.csv('data/gisaid_terminal_4-fold.csv')  # 0-based, true 4_fold

# ss_annot = read.csv('data/secondary_structure_on_genome.csv')  # 1-based
pre_genome_annot <- read.csv('data/gemone_structures_annot_with_UTRs.csv')  # 1-based

pre_genome_annot$IsStem[is.na(pre_genome_annot$IsStem)] <- -1  # fill NA
genome_annot <- pre_genome_annot[c('Pos', 'IsStem', 'IsTRS')]

# prepare gisaid
gis$Pos = gis$Pos + 1  # 1-based
names(gis)[7:8] = c('RefNuc', 'AltNuc')
gis = subset(gis, RefNuc != '-' & AltNuc != '-')
gis$AaSub = ifelse(gis$child_aa == gis$parent_aa, 'S', 'NS')
head(gis)

# only leaves of tree to catch mutagenesis | already done for 4-fold
# gis = subset(gis, ! startsWith(gis$child_node, '#'))

# -------------------------
# MutSpec for all positions
draw_mutspec(df = gis, ann = ann.ff, n_mutspec = 12, 
             mut_type = 'S',
             spec_type = 'general')

# MutSpec for stem positions
# gis.merged <- merge(gis, genome_annot, by='Pos')
gis.merged <- merge(gis, genome_annot[c(1, 3)], by='Pos')
gis.stems = subset(gis.merged, IsStem == 1)
gis.free = subset(gis.merged, IsStem == 0)

ann.merged = merge(ann.ff, genome_annot, by='Pos')
ann.stems = subset(ann.merged, IsStem == 1)
ann.free = subset(ann.merged, IsStem == 0)

draw_mutspec(df = gis.stems, ann = ann.stems, n_mutspec = 12, 
             mut_type = 'S',
             spec_type = 'general', 
             label = 'Stem')

draw_mutspec(df = gis.free, ann = ann.free, n_mutspec = 12, 
             mut_type = 'S',
             spec_type = 'general', 
             label = 'Free')

# MutSpec for TRS and RS(?)
# gis.merged <- merge(gis, genome_annot, by.x='pos', by.y='Pos')
# ann.merged = merge(ann, genome_annot, by='Pos')

gis.trs <- subset(gis.merged, IsTRS == 1)
ann.trs <- subset(ann.merged, IsTRS == 1)

draw_mutspec(df = gis.trs, ann = ann.trs, n_mutspec = 12, 
             mut_type = 'S',
             spec_type = 'general', 
             label = 'TRS')

# MutSpec for genes in stems
draw_mutspec(df = gis.stems, ann = ann.stems, n_mutspec = 12, 
             mut_type = 'S',
             spec_type = 'genes', 
             label = 'Stem')

# MutSpec for genes in all pos
draw_mutspec(df = gis.free, ann = ann.free, n_mutspec = 12, 
             mut_type = 'S',
             spec_type = 'genes',
             label = 'Free')


# sort(table(subset(pre_genome_annot, IsTRS == 1)$RefNuc))
# sort(table(gis.trs$GenName), decreasing = TRUE)
# 
# sort(table(
#   subset(gis.trs, RefNuc == 'T' & AltNuc == 'C')$pos), 
#   decreasing = TRUE
# )
# 
# sort(table(
#   subset(gis.trs, pos == 27385)$RefNuc),
#   decreasing = TRUE
# )




#This function takes at the entres two data frames, df - data frame with information from nextstrain, ann - table with mutations, neighbours, aminoacids and etc.
# Possible variable values for drawing n_mutspec = 12/192, mut_type = 'S'/'NS', spec_type = 'general'/'genes'/'data', region_type = 'translated'/'untranslated'
# draw_mutspec(df = gis.sample, ann = ann.sample, n_mutspec = 12, 
#              mut_type = 'S',
#              spec_type = 'general', 
# )


read.csv('data/MutSpecTables/07.MutSpec12_ForFullGenome.csv')
read.csv('data/MutSpecTables/07.MutSpec12_ForStem.csv')


