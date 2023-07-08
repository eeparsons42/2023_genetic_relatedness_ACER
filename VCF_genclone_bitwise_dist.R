library(vcfR)
library(poppr)
library(adegenet)
library(ape)
library(ggplot2)
library(knitr)
library(RColorBrewer)



vcf<-read.vcfR("C:/Users/parsonsee/Desktop/Genetic analysis practice/dup_removed_DB_snps.vcf.gz")

gt <- extract.gt(vcf, element = "GT", as.numeric=FALSE)
myMiss <- apply(gt, MARGIN = 2, function(x){ sum(is.na(x)) })
myMiss <- (myMiss/nrow(vcf))*100
miss <- data.frame(myMiss)
write.table(miss,file="missingData.txt",quote =FALSE)

#Convert the VCF file into a format compatible with the Poppr package
#genind object
gl <- vcfR2genlight(vcf, n.cores=2)


gind <- new("genind",(as.matrix(gl))) #all SNPs
#add population information to the genind object
poptab<-read.table("popInfo_p1.txt",check.names=FALSE, header=T, na.strings = c("", "NA"))
gind@pop <- as.factor(poptab$region)
gind

#convert genind to genclone object
gclo<-as.genclone(gind)
gclo

#calculate the bitwise distance between individuals
xdis<-bitwise.dist(gclo, missing_match=FALSE) 
#similar to Provesti's distance

##Create a phylogeny of samples based on distance matrices
#colors
cols <- palette(brewer.pal(n=12, name = 'Set3'))
set.seed(999)
theTree <- gclo %>%
  aboot(dist = provesti.dist, sample = 10, tree = "nj",
        cutoff = 50, quiet = TRUE) %>%
  ladderize() # organize branches by clade

plot.phylo(theTree, tip.color = cols[gclo$pop],label.offset = 0.0125,cex=0.7,
           font=2, lwd=4,align.tip.label= F,no.margin = T)
add.scale.bar(0,0.95,length = 0.05, cex=0.65,lwd=3) # add a scale bar showing 5% difference.
nodelabels(theTree$node.label, cex=.5, adj = c(1.5, -0.1), frame = "n", font=3,
           xpd = TRUE)
#substrees of the full tree
l<-subtrees(theTree)
#A. palmata clones
plot(l[[113]])
add.scale.bar(0,0.95,length = 0.01, cex=0.65,lwd=3)
