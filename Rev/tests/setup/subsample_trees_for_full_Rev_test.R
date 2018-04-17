setwd("~/git_repos/jfaulkner_phylocode/Rev/tests/")

# Thin bottleneck
tree <- read.tree("tree_only/data/bottleneck.tre")
taxa <- read.table("tree_only/data/bottlneck.taxa.txt",sep=" ",header=TRUE,stringsAsFactors=FALSE)

keep <- round(seq(1,dim(taxa)[1],length.out=100))
tree <- drop.tip(tree,taxa$taxon[-keep])

write.tree(tree,"full_model/data/bottleneck.tre",digits=30)
write.table(taxa[keep,],"full_model/data/bottleneck.taxa.txt",sep=" ",quote=FALSE,row.names=FALSE)

# Thin mex hat
tree <- read.tree("tree_only/data/mex_hat.tre")
taxa <- read.table("tree_only/data/mex_hat.taxa.txt",sep=" ",header=TRUE,stringsAsFactors=FALSE)

keep <- round(seq(1,dim(taxa)[1],length.out=100))
tree <- drop.tip(tree,taxa$taxon[-keep])

write.tree(tree,"full_model/data/mex_hat.tre",digits=30)
write.table(taxa[keep,],"full_model/data/mex_hat.taxa.txt",sep=" ",quote=FALSE,row.names=FALSE)

# Thin bottleneck
tree <- read.tree("tree_only/data/broken_exponential.tre")
taxa <- read.table("tree_only/data/broken_exponential.taxa.txt",sep=" ",header=TRUE,stringsAsFactors=FALSE)

keep <- round(seq(1,dim(taxa)[1],length.out=100))
tree <- drop.tip(tree,taxa$taxon[-keep])

write.tree(tree,"full_model/data/broken_exponential.tre",digits=30)
write.table(taxa[keep,],"full_model/data/broken_exponential.taxa.txt",sep=" ",quote=FALSE,row.names=FALSE)

# Thin bottleneck
tree <- read.tree("tree_only/data/GP.tre")
taxa <- read.table("tree_only/data/GP.taxa.txt",sep=" ",header=TRUE,stringsAsFactors=FALSE)

keep <- round(seq(1,dim(taxa)[1],length.out=100))
tree <- drop.tip(tree,taxa$taxon[-keep])

write.tree(tree,"full_model/data/GP.tre",digits=30)
write.table(taxa[keep,],"full_model/data/GP.taxa.txt",sep=" ",quote=FALSE,row.names=FALSE)
