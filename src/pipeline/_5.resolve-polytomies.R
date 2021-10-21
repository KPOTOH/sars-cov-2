#!/mnt/lustre/genkvg/anaconda2/bin/Rscript
#USAGE: ./5.resolve-polytomies.R fastaMultipleAlignment.tre-simple binarySimplifiedTree

# http://blog.phytools.org/2017/06/generating-set-of-random-resolutions-of.html

library(phytools)
library(phangorn)
library(ape)
args = commandArgs(trailingOnly=TRUE)

resolveRandom<-function(tree){
while(!is.binary(tree)){
    nodes<-1:tree$Nnode+Ntip(tree)
    Nchildren<-function(node,tree) length(Children(tree,node))
    nchilds<-sapply(nodes,Nchildren,tree=tree)
    node<-nodes[which(nchilds>2)[1]]
    tree<-sample(resolveNode(tree,node),1)[[1]]}
tree
}

tree<-read.tree(args[1])
treeResolved<-resolveRandom(tree)
write.tree(treeResolved, args[2])
