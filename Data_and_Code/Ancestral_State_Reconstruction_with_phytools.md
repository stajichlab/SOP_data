# Ancestral State Reconstruction (ASR) using [phytools](https://github.com/liamrevell/phytools)
## 1. Identify the type of data you are interested in analyzing
- Continuous characters (characters on a continuous scale)
- Discrete characters (qualitative features or characteristics we can count)

## Tree Reconstruction

## Ancestral State Reconstruction

### Discrete characters

#### Example of discrete character script to create two PDFs, each for a tree with 1 trait

```R
## Conidiobolus ASR - two traits, four conditions/trait
## 09/30/21
## Kelsey Aadland
## Script to create two plots with ASR values for each trait
  # both are scaled to be exported as PDFs at width=10" height=15"
  
library(phytools)

#read in the tree and extant state data
tree <- read.tree("NieGenomeComb.prune_root.tre")
tree <- force.ultrametric(tree, method=c("nnls","extend"))
x <- as.matrix(read.csv("Conidiobolus_ASR.prune.csv",row.names=1))[,1]
y <- as.matrix(read.csv("Conidiobolus_ASR.prune.csv",row.names=1))[,2]
xz <- factor(x=as.character(x),levels=c("0","1","2","3"))
yz <- factor(x=as.character(y),levels=c("0","1","2","3"))

#set pdf size for trait 1
pdf(file="Conidiobolus_ASR_trait1.pdf",width=10,height=15)
#estimate ancestral states for trait 1
objx<-rerootingMethod(tree,x,model="ER")
# plot tree
plot(objx,fsize=0.5,cex=c(0.5,0.5))
dev.off()

#set pdf size for trait 2
pdf(file="Conidiobolus_ASR_trait2.pdf",width=10,height=15)
#estimate ancestral states for trait 2
objy<-rerootingMethod(tree,y,model="ER")
# plot tree
plot(objy,fsize=0.5,cex=c(0.5,0.5))
dev.off()
```
*For more tree visualization methods, check out* [Exercise 15: Plotting methods for phylogenies & comparative data in R](http://www.phytools.org/Cordoba2017/ex/15/Plotting-methods.html)

## Citations

* Revell, L.J. (2012), phytools: an R package for phylogenetic comparative biology (and other things). Methods in Ecology and Evolution, 3: 217-223. https://doi.org/10.1111/j.2041-210X.2011.00169.x

## Resources

- Phytools [documentation](https://cran.r-project.org/web/packages/phytools/phytools.pdf) (Updated 09/01/2022)
- Phytools [blog](http://blog.phytools.org/)
  - Phytools blog entry - [Plotting multiple pie charts at nodes](http://blog.phytools.org/2016/10/plotting-multiple-pie-charts-at-nodes.html)
- "Exercise 15: [Plotting methods for phylogenies & comparative data in R"](http://www.phytools.org/Cordoba2017/ex/15/Plotting-methods.html)
- 
