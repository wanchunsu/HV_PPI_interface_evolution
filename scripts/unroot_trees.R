# install required package if not yet installed
if (! "ape" %in% installed.packages()){
  install.packages(i, dependencies = TRUE)
}
require("ape")

# Load tree
tree_fi = "../data/binned_dnds/binned_dnds.species_tree.nwk"
tr <- ape::read.tree(tree_fi)

# Unroot tree
unrooted_tr <- ape::unroot(tr)

# Write unrooted tree to new file
unrooted_tree_fi = "../data/binned_dnds/unrooted_binned_dnds.species_tree.nwk"
write.tree(unrooted_tr, unrooted_tree_fi)