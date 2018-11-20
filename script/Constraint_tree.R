library("ape")

#contraint_tree <- read.tree(text="(((((((((Lycium_villosum,Lycium_hirsutum),Lycium_bosciifolium),Lycium_shawii),Lycium_europaeum),(((Lycium_barbarum,Lycium_truncatum),(Lycium_dasystemum,Lycium_ruthenicum)),(Lycium_chinense,Lycium_yunnanense))),Lycium_pumilum),((Lycium_tenue,Lycium_ferocissimum),Lycium_oxycarpum)),((Lycium_chilense,Lycium_carolinianum),Lycium_andersonii)),Nolana_werdermannii);")

contraint_tree <- read.tree("../results/Lycium_constraint_tree-Miller_et_al_2011.tre")

pdf("../results/Lycium_constraint_tree-Miller_et_al_2011.pdf")
plot(contraint_tree)
dev.off()

write.tree(contraint_tree, "../results/Lycium_constraint_tree-Miller_et_al_2011.tre")