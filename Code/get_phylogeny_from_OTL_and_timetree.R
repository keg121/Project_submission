#!/usr/bin/env Rscript

library(ape)
library(geiger)
library(rotl)

## load in data
df <- read.csv("TPC_parameter_estimates.csv", row.names=1)

## remove unmatched species/ species not found in tree of life topology

# no correct synonym
df <- df[!grepl("Rhynchaenus flagellum", df$Species),] 
df <- df[!grepl("Encoptolophus sordidus costalis", df$Species),]
df <- df[!grepl("Pterohelaeus spp.", df$Species),]
df <- df[!grepl("Cerotalis spp.", df$Species),]
df <- df[!grepl("Carenum spp.", df$Species),]
df <- df[!grepl("Glossina morsitans orientalis", df$Species),]

#not in OTL 
df <- df[!grepl("Palirhoeus eatoni", df$Species),]
df <- df[!grepl("Hophlosphyrum griseus", df$Species),]
df <- df[!grepl("Helius waitei", df$Species),]
df <- df[!grepl("Scarabaeus gariepinus", df$Species),]
df <- df[!grepl("Scarabaeus hippocrates", df$Species),]



# The unique species names from your dataset should be in this 'species_names' vector.
species_names <- unique(df[[1]])

# This function matches the species names to those in the Open Tree of Life.
# The function returns a data frame (resolved_names). You have to manually 
# examine it to see if the names it found are indeed synonyms or not. 
#
# If a name is wrongly matched, then you may want to remove that species 
# from your original vector (species_names).
resolved_names <- tnrs_match_names(species_names)


# This returns a tree topology from the Open Tree of Life based on the 
# species in the resolved_names data frame.
#
# It may return an error if one of the species is not present in their 
# topology tree. If that is the case, you have to manually remove that 
# species from your original vector (species_names).
tree <- tol_induced_subtree(ott_ids=ott_id(resolved_names))

# The Open Tree of Life topology always has an ID number at the end of 
# each species' name. This line below removes that ID.
tree$tip.label <- gsub('_ott\\d+$', '', tree$tip.label)

# Then, you need to go to timetree.org and upload your list of species 
# names (one per line). This will give you a time-calibrated tree of 
# as many species as possible.

# As before, it may replace some of the species with synonyms. You need 
# to manually examine the replacements to make sure that the synonyms 
# are valid. If any of them is not, it is best to remove the species 
# that triggered a wrong replacement from the list you upload to timetree.org.
#
# Do not worry if some of the species in your list are not present in 
# this tree. Once you are happy with the tree that timetree.org shows, 
# download it and read it in R with the command below.
timetree_tree <- read.tree("species.nwk")

# Finally, the function below will transfer ages from the tree you got 
# from timetree.org to the OTL topology.
#
# Any nodes that are incompatible between the two trees or for which you 
# do not have any age information will be scaled using the treePL program. 
# You need to have that program manually installed and in your PATH 
# (i.e., if you  type treePL in the terminal, the program should run).
final_calibrated_timetree <- congruify.phylo(
	timetree_tree, tree, scale = 'treePL'
)

# Finally, this function writes the calibrated OTL tree to a file.
write.tree(final_calibrated_timetree$phy, file = "final_calibrated_tree.phy")
