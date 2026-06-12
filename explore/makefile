tree_original = "../config/polished_tree.nwk"
tree_filter   = "../results/core_trees/polished_from_filtering.nwk"
tree_gubbins  = "../results/gubbins/gubbins.final_tree.tre"

out_root = figs/n04_treeknit

auspice_view: original_vs_filter original_vs_gubbins filter_vs_gubbins

original_vs_filter:
	mkdir -p $(out_root)/original_vs_filter
	treeknit $(tree_original) $(tree_filter)  -o $(out_root)/original_vs_filter/treeknit  --auspice-view

original_vs_gubbins:
	mkdir -p $(out_root)/original_vs_gubbins
	treeknit $(tree_original) $(tree_gubbins) -o $(out_root)/original_vs_gubbins/treeknit --auspice-view

filter_vs_gubbins:
	mkdir -p $(out_root)/filter_vs_gubbins
	treeknit $(tree_filter)   $(tree_gubbins) -o $(out_root)/filter_vs_gubbins/treeknit  --auspice-view
