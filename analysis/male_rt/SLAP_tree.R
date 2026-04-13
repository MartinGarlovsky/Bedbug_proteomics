library(tidyverse)
library(ape)        # for phylogenetic analysis
library(phangorn)   # for tree building methods
library(seqinr)     # for sequence handling
library(ggtree)     # optional, for publication-quality trees

# Read the MAFFT alignment file
alignment_file <- "output/male_rt/SLAP_stuff/M17SLAP3_SLAP12_GrSm_aligned.fasta"
alignment <- read.FASTA(alignment_file, type = "AA")

# Convert to phyDat object (required for phangorn)
alignment_phyDat <- as.phyDat(alignment)#, type = "AA")

# Calculate distance matrix using "JTT" (Jones-Taylor-Thornton)
dist_matrix <- dist.ml(alignment_phyDat, model = "JTT")

# Start with NJ tree as initial tree
init_tree <- nj(dist_matrix)

# Root the tree to LAP3 outgroup
init_tree <- root(init_tree, outgroup = "P00727", resolve.root = TRUE)

ggtree(init_tree) + geom_tiplab()

# Optimize tree using Maximum Likelihood
ml_tree <- pml(init_tree, alignment_phyDat, model = "JTT")

# # Optimize tree topology and parameters # this takes time to run
# ml_tree_optimized <- optim.pml(ml_tree,
#                                model = "JTT",
#                                optNni = TRUE,      # optimize topology
#                                optBf = TRUE,       # optimize base frequencies
#                                optQ = TRUE,        # optimize rate matrix
#                                optGamma = TRUE,    # optimize gamma parameter
#                                rearrangement = "NNI")  # tree rearrangement method

#write_rds(ml_tree_optimized, "output/male_rt/SLAP_stuff/mafft_output/ml_tree_optimized_LAP3.rds")
ml_tree_optimized <- read_rds("output/male_rt/SLAP_stuff/mafft_output/ml_tree_optimized_LAP3.rds")

# Extract the best tree
best_tree <- ml_tree_optimized$tree

# Bootstrap analysis (assess branch support)
# bs_trees <- bootstrap.pml(ml_tree_optimized,
#                           bs = 1000,
#                           optNni = TRUE)

#write_rds(bs_trees, "output/male_rt/SLAP_stuff/mafft_output/bs_trees_LAP3.rds")
bs_trees <- read_rds("output/male_rt/SLAP_stuff/mafft_output/bs_trees_LAP3.rds")

# Add bootstrap values to tree
best_tree_bs <- plotBS(best_tree, bs_trees, type = "unrooted",
                       bs.col = "red", bs.adj = c(1.2, 0.5))

#write_rds(best_tree_bs, "output/male_rt/SLAP_stuff/mafft_output/best_tree_bs.rds")
