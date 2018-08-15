## This script is for subsetting the simulated data for the full_model analyses

#setwd("~/git_repos/jfaulkner_phylocode/Rev/tests/")

# Thin bottleneck
tree <- read.tree("tree_only/data/bottleneck.tre")
taxa <- read.table("tree_only/data/bottleneck.taxa.txt",header=TRUE,stringsAsFactors=FALSE)

by_increasing_age <- order(taxa[,2])
taxa <- taxa[by_increasing_age,]

keep <- round(seq(1,dim(taxa)[1],length.out=100))
tree <- drop.tip(tree,taxa$taxon[-keep])

write.tree(tree,"full_model/data/bottleneck.tre",digits=30)
write.table(taxa[keep,],"full_model/data/bottleneck.taxa.txt",quote=FALSE,row.names=FALSE)

# Thin mex hat
tree <- read.tree("tree_only/data/mex_hat.tre")
taxa <- read.table("tree_only/data/mex_hat.taxa.txt",header=TRUE,stringsAsFactors=FALSE)

by_increasing_age <- order(taxa[,2])
taxa <- taxa[by_increasing_age,]

keep <- round(seq(1,dim(taxa)[1],length.out=100))
tree <- drop.tip(tree,taxa$taxon[-keep])

write.tree(tree,"full_model/data/mex_hat.tre",digits=30)
write.table(taxa[keep,],"full_model/data/mex_hat.taxa.txt",quote=FALSE,row.names=FALSE)

# Thin Broken Exponential
tree <- read.tree("tree_only/data/broken_exponential.tre")
taxa <- read.table("tree_only/data/broken_exponential.taxa.txt",header=TRUE,stringsAsFactors=FALSE)

by_increasing_age <- order(taxa[,2])
taxa <- taxa[by_increasing_age,]

keep <- round(seq(1,dim(taxa)[1],length.out=100))
tree <- drop.tip(tree,taxa$taxon[-keep])

write.tree(tree,"full_model/data/broken_exponential.tre",digits=30)
write.table(taxa[keep,],"full_model/data/broken_exponential.taxa.txt",quote=FALSE,row.names=FALSE)

# Thin GP
tree <- read.tree("tree_only/data/GP.tre")
taxa <- read.table("tree_only/data/GP.taxa.txt",header=TRUE,stringsAsFactors=FALSE)

by_increasing_age <- order(taxa[,2])
taxa <- taxa[by_increasing_age,]

keep <- round(seq(1,dim(taxa)[1],length.out=100))
tree <- drop.tip(tree,taxa$taxon[-keep])

write.tree(tree,"full_model/data/GP.tre",digits=30)
write.table(taxa[keep,],"full_model/data/GP.taxa.txt",quote=FALSE,row.names=FALSE)
