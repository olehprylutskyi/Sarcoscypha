#reset ctrl+shift+f10
rm(list=ls())
library(ape)
library(ips)
sarco.dna<-read.dna(file="sarco_seq.fasta", format = "fasta")
sarco.dna
class(sarco.dna)
sarco.mafft <- mafft(x=sarco.dna, method="localpair", maxiterate=100, options="--adjustdirection", exec="C:/Users/Mycolab-2017/Desktop/mafft-win/mafft")
sarco.mafft
sarco.mafft.ng<-deleteGaps(x=sarco.mafft,gap.max=nrow(sarco.mafft)/4)
sarco.mafft.ng
class(sarco.mafft.ng)
image.DNAbin(sarco.mafft.ng)
library(ape) # Analysis of phylogenetics and evolution
library(hierfstat) # Hierarchical F-statistics
library(corrplot) # Visualization of correlation matrix
# Create the distance matrix
sarco.dist <- dist.dna(x=sarco.mafft.ng, model="TN93")
# Check the resulting distance matrix
sarco.dist
class(sarco.dist)
dim(as.matrix(sarco.dist))
# Create another distance matrix
sarco.dist1 <- dist.dna(x=sarco.mafft.ng, model="F81")
# Check it
sarco.dist1
class(sarco.dist1)
dim(as.matrix(sarco.dist1))
class(sarco.dist1)
table.paint(df=as.data.frame(as.matrix(sarco.dist1)), cleg=0, clabel.row=0.5, clabel.col=0.5)
# Same visualization, colored
# heatmap() reorders values
# because by default it plots
# also dendrograms on the edges
heatmap(x=as.matrix(sarco.dist1), Rowv=NA, Colv=NA, symm=TRUE)
# According to distance used
# How to use hierarchical clustering
#This is very basic function to make dendrogram
?hclust
plot(hclust(d=sarco.dist1, method="complete"))
# Calculate it
# Saving as phylo object (and not hclust) gives more
# possibilities for further plotting and manipulations
sarco.upgma <- as.phylo(hclust(d=sarco.dist1, method="average"))
plot.phylo(x=sarco.upgma, cex=0.75)
title("UPGMA tree")
# Test quality - tests correlation of original distance in the matrix
# and reconstructed distance from hclust object
plot(x=as.vector(sarco.dist1), y=as.vector(as.dist(
cophenetic(sarco.upgma))), xlab="Original pairwise distances",
ylab="Pairwise distances on the tree", main="Is UPGMA appropriate?", pch=20, col=transp(col="black",
                                     alpha=0.1), cex=2)
# Add correlation line
abline(lm(as.vector(as.dist(cophenetic(sarco.upgma)))~
               as.vector(sarco.dist1)), col="red")
#Neighbor-Joining tree
# Test tree quality - plot original vs. reconstructed distance
sarco.nj <- nj(sarco.dist1)
class(sarco.nj)
plot(as.vector(sarco.dist1), as.vector(as.dist(cophenetic(sarco.nj))),
xlab="Original distance", ylab="Reconstructed distance")
abline(lm(as.vector(sarco.dist1) ~
as.vector(as.dist(cophenetic(sarco.nj)))), col="red")
cor.test(x=as.vector(sarco.dist1), y=as.vector(as.dist(cophenetic
(sarco.nj))), alternative="two.sided") # Testing the correlation
# Linear model for above graph
summary(lm(as.vector(sarco.dist1) ~
               as.vector(as.dist(cophenetic(sarco.nj))))) # Prints summary text
# Plot a basic tree - see ?plot.phylo for details
plot.phylo(x=sarco.nj, type="phylogram")
plot.phylo(x=sarco.nj, type="cladogram", edge.width=2)
plot.phylo(x=sarco.nj, type="fan", edge.width=2, edge.lty=2)
plot.phylo(x=sarco.nj, type="radial", edge.color="red", edge.width=2,
edge.lty=3, cex=2) 
# tree based on DNA sequences
sarco.tree <- nj(X=sarco.dist1)
plot.phylo(x=sarco.tree, type="unrooted", show.tip=FALSE)
title("Unrooted NJ tree")
# Read annotations
sarco.annot <- read.csv("SarcoData.csv", header=TRUE, row.names=1)
head(sarco.annot) # See result
# Colored tips
sarco.pal <- colorRampPalette(topo.colors(length(levels(as.factor(sarco.annot[["species"]])))))
# Tip labels
tiplabels(text=sarco.annot$accession, bg=num2col(sarco.annot$accession,
                                           col.pal=sarco.pal), cex=0.75)
# Calculate bootstrap
sarco.boot <- boot.phylo(phy=sarco.tree.rooted, x=sarco.dna, FUN=function(EEE) root.phylo(nj(dist.dna(EEE, model="TN93")), outgroup=1),B=1000)
#Rooting and unrooting trees
plot.phylo(sarco.nj)
print.phylo(sarco.nj)
sarco.nj.rooted <- root.phylo(phy=sarco.nj, interactive=TRUE)
plot.phylo(sarco.nj.rooted)
#Maximum parsimony
library(phangorn)
# Conversion to phyDat for phangorn
sarco.phydat <-as.phyDat(sarco.mafft.ng)
# Prepare starting tree
sarco.tre.ini <- nj(dist.dna(x=sarco.mafft.ng, model="raw"))
?parsimony # Parsimony details
# Returns maximum parsimony score
parsimony(tree=sarco.tre.ini,data=sarco.phydat)
# Optimisation - returns MP tree
sarco.tre.pars <- optim.parsimony(tree=sarco.tre.ini,data=sarco.phydat)
plot.phylo(x=sarco.tre.pars,type="cladogram", edge.width=2)
