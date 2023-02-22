
library(phytools)
library(phangorn)

x = 'cat'
y = 'squirrel'

z = paste(x, y, sep = '')
z


setwd('~/OneDrive - Imperial College London/TreeBuilding/230203filter/dating/')

getwd()

tre <- read.nexus("aa_crown_230207.timetree.nex")

# Check if tree is rooted: BAMM tools needs rooted tree
is.rooted(tre)
# If it needs to be rooted, check node number to root it from
plot.phylo(tre)
nodelabels()
tiplabels()
# Root tree
tre <- root(tre, 1)


# Make binary
is.binary(tre)            # check if it's fully resolved
tre <- multi2di(tre)

# Remove 0 branch lengths
any(tre$edge.length == 0) # check if it has any 0 branch lengths
tre <- dispRity::remove.zero.brlen(tre)

# Make ultrametric (already should be, but might be slightly off)
is.ultrametric(tre)
N <- Ntip(tre)
root_node <- N + 1
root_to_tip <- dist.nodes(tre)[1:N, root_node]
age_difference <- max(root_to_tip) - root_to_tip
tre_extend <- tre
tip_edges <- tre_extend$edge[, 2] <= Ntip(tre_extend)
tre_extend$edge.length[tip_edges] <- tre_extend$edge.length[tip_edges] + age_difference
is.ultrametric(tre_extend)
write.tree(tre_extend, file="ultrametric_aa.tre")

tre <- read.tree('ultrametric_aa.tre')

