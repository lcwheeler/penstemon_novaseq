#Script used to transform gene trees to ultrametric.
#These will be used to visualize discordance in DensiTree

library(ape)

#read in trees
setwd('~/Desktop/')
intrees <- read.tree("combined_10kbwindowtrees.tre")


#root trees
x <- root(intrees, outgroup = "mon_61-7_S440", resolve.root = TRUE)


#make trees ultrametric
for (i in 1:length(x)){
  x[[i]] <- chronos(x[[i]], lambda = 0)
}


#write trees to outfile
write.tree(x, file = "ultrametric_combined_10kbwindowtrees.tre")
