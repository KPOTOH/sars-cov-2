# install.packages("ggplot2")
# http://www.sthda.com/english/wiki/ggplot2-barplots-quick-start-guide-r-software-and-data-visualization

library('ggplot2')
library(LaplacesDemon)

setwd("~/sars-cov-2/")
rm(list=ls(all=TRUE))

AddNormedFr <- function(df){
  df$OverPres <- df$ExpFr / min(df$ExpFr)
  df$NormedObsFr <- df$ObsFr / df$OverPres
  ura_names <- c("A>C", "A>G", "A>U", "C>A", "C>G", "C>U", 
                 "G>A", "G>C", "G>U", "U>A", "U>C", "U>G")
  df$NucSubst <- ura_names
  return(df)
}

bootstraping = function(df){
  return(df[sample(nrow(df), replace = TRUE),])
}


nxt_mutspec <- read.csv('data/share/MutSpecTables/07.MutSpec12_Nextstrain_FF.csv')
full_mutspec <- read.csv('data/MutSpecTables/07.MutSpec12_ForAllPos.csv')
stem_mutspec <- read.csv('data/MutSpecTables/07.MutSpec12_ForStem.csv')
free_mutspec <- read.csv('data/MutSpecTables/07.MutSpec12_ForFree.csv')
trs_mutspec <- read.csv('data/MutSpecTables/07.MutSpec12_ForTRS_RS.csv')

nxt_mutspec <- AddNormedFr(nxt_mutspec)
full_mutspec <- AddNormedFr(full_mutspec)
free_mutspec <- AddNormedFr(free_mutspec)
stem_mutspec <- AddNormedFr(stem_mutspec)
trs_mutspec <- AddNormedFr(trs_mutspec)

KLD(full_mutspec$MutSpec, nxt_mutspec$MutSpec)
KLD(stem_mutspec$MutSpec, free_mutspec$MutSpec)

#################################################
# Plot nextstrain versus gisaid 
#################################################
full_mutspec$Набор_замещений = 'GISAID'
nxt_mutspec$Набор_замещений = 'Nextstrain'
  
df = rbind(full_mutspec, nxt_mutspec)

p <- ggplot(data=df, aes(x=NucSubst, y=MutSpec, fill=Набор_замещений)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()

# Change color by groups
# Add error bars
p + labs(x="Замещения", 
         y = "Мутационный спектр") +
  scale_fill_manual(values=c('white', 'gray49'))+
  theme_classic()

M <- as.table(rbind(nxt_mutspec$NormedObsFr, full_mutspec$NormedObsFr))
M <- round(M)
(Xsq <- chisq.test(M, simulate.p.value = TRUE))  # Prints test summary

#################################################
# Plot gisaid subsamples according to ss
#################################################

full_mutspec$Позиции = 'Все'
stem_mutspec$Позиции = 'Спаренные'
free_mutspec$Позиции = 'Неспаренные'
trs_mutspec$Позиции = 'ТRS & RS'

# Chi-square test
M <- as.table(rbind(stem_mutspec$NormedObsFr, free_mutspec$NormedObsFr))
M <- round(M)
(Xsq <- chisq.test(M, simulate.p.value = TRUE))

# Fisher test
N <- length(stem_mutspec$NormedObsFr)
for (i in 1:N) {
  stemCnt <- stem_mutspec$NormedObsFr[i]
  freeCnt <- free_mutspec$NormedObsFr[i]
  
  others_stemCnt <- sum(stem_mutspec$NormedObsFr[1:N != i])
  others_freeCnt <- sum(free_mutspec$NormedObsFr[1:N != i])
  
  fdf <- round(matrix(c(stemCnt, others_stemCnt, freeCnt, others_freeCnt), nrow = 2))
  fresult <- fisher.test(fdf)
  print(fresult$p.value)
}

df = rbind(full_mutspec, stem_mutspec, free_mutspec)
df = rbind(stem_mutspec, trs_mutspec)

# Change the colors manually
p <- ggplot(data=df, aes(x=NucSubst, y=MutSpec, fill=Позиции)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()

# Change color by groups
# Add error bars
p + labs(x="Замещения", 
         y = "Мутационный спектр") +
  # scale_fill_manual(values=c('white','lightgray', 'gray49'))+
  scale_fill_manual(values=c('lightgray', 'black'))+
  theme_classic()

# ------------------------------

genes <- c('ORF1ab', 'M')
for (gname in genes){
  stem_mutspec.genes <- read.csv(sprintf('data/MutSpecTables/07.MutSpec12_%s_%s.csv', 'Stem', gname))
  full_mutspec.genes <- read.csv(sprintf('data/MutSpecTables/07.MutSpec12_%s_%s.csv', 'Free', gname))
  
  stem_mutspec.genes$Type <- 'Stem'
  full_mutspec.genes$Type <- 'Free'
  
  both_mutspecs = rbind(stem_mutspec.genes, full_mutspec.genes)
  
  p <- ggplot(data=both_mutspecs, aes(x=NucSubst, y=MutSpec, fill=Type)) +
    geom_bar(stat="identity", color="black", position=position_dodge())+
    theme_minimal()
  
  p <- p + labs(title=sprintf("MutSpec for gene %s", gname), x="Substitutions", y = "MutSpec") +
    scale_fill_manual(values=c('white','gray40'))+
    theme_classic()
  show(p)
}

# Full gemone without TRS bars
df = rbind(full_mutspec, stem_mutspec)

# Change the colors manually
p <- ggplot(data=df, aes(x=NucSubst, y=MutSpec, fill=Type)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()

# Change color by groups
# Add error bars
p + labs(title="MutSpec on subsamples of substitutions", x="Substitutions", y = "MutSpec") +
  scale_fill_manual(values=c('white','lightgray'))+
  theme_classic()


