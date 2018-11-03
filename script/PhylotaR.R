#install.packages('phylotaR')
library("phylotaR")
library("ggplot2")
source("./tools/selection_tools.R")#scripts from :https://github.com/AntonelliLab/phylotaR_demo
# dir.create('Lycium')
# wd <- getwd()
# blast_dir="/apps/blast/ncbi/2.6.0/bin/"
# setup(wd = wd, txid = 24646, ncbi_dr = blast_dir)
# run(wd = wd)

wd <- "../data/PhylotaR_test"
phylota <- read_phylota(wd)

# drop all but first ten
Lycium <- drop_clstrs(phylota, phylota@cids[1:10])

#Lycium2 <- drop_clstrs(phylota, phylota@cids[1:7])
# get genus-level taxonomic names
species_txids <- get_txids(Lycium, txids = Lycium@txids, rnk = 'species')

species_txids <- unique(species_txids)

# dropping missing
species_txids <- species_txids[species_txids !=  '']
species_nms <- get_tx_slot(Lycium, species_txids, slt_nm = 'scnm')

# make alphabetical for plotting
species_nms <- sort(species_nms, decreasing = TRUE)
# generate geom_object
p <- plot_phylota_pa(phylota = Lycium, cids = Lycium@cids, txids = species_txids,
                     txnms = species_nms)
# plot
print(p)

q <- plot_phylota_treemap(phylota = Lycium, txids = species_txids, txnms = species_nms,
                          area = 'nsq', fill = 'ncl')
print(q)


cid <- Lycium@cids[[1]]
# extract its sequence IDs from Phylota object
sids <- Lycium[[cid]]@sids
# design sequence names for definition line
sid_txids <- get_sq_slot(phylota = Lycium, sid = sids, slt_nm = 'txid')
sid_spnms <- get_tx_slot(phylota = Lycium, txid = sid_txids, slt_nm = 'scnm')
sq_nms <- paste0(sid_spnms, ' | ', sids)
# write out
write_sqs(phylota = Lycium, "../data/Lycium.fasta", sq_nm = sq_nms, sid = sids)

##############################################################################
all_cls <- phylota
# RECREATE PhyLoTa TABLE
# get phylogenetically inform. clusters
n_taxa <- get_clstr_slot(all_cls, cid = all_cls@cids, slt_nm = 'ntx')
cltyps <- get_clstr_slot(all_cls, all_cls@cids, slt_nm = 'typ')
#picls <- cltyps %in% c('subtree', 'direct') & n_taxa >= 4
picls <- n_taxa >= 4
picls <- all_cls@cids[picls]
# build table
prnts <- get_clstr_slot(all_cls, cid = picls, slt_nm = 'prnt')
n_taxa <- get_clstr_slot(all_cls, cid = picls, slt_nm = 'ntx')
n_sqs <- get_clstr_slot(all_cls, cid = picls, slt_nm = 'nsqs')
# n_gnra <- sapply(picls, function(x) {
#   length(unique(get_txids(all_cls, cid = x, rnk = 'genus')))
# })

n_species <- sapply(picls, function(x) {
  length(unique(get_txids(all_cls, cid = x, rnk = 'species')))
})

#max seq length
mxsqlngs <- sapply(picls, function(x) {
  max(get_sq_slot(all_cls, cid = x, slt_nm = 'nncltds'))
})
#min seq length
mnsqlngs <- sapply(picls, function(x) {
  min(get_sq_slot(all_cls, cid = x, slt_nm = 'nncltds'))
})


lngst <- sapply(picls, function(x) {
  sids <- all_cls[[x]]@sids
  lngths <- get_sq_slot(all_cls, sid = sids, slt_nm = 'nncltds')
  sids[which.max(lngths)]
})

deflns <- get_sq_slot(all_cls, sid = lngst, slt_nm = 'dfln')
mads <- sapply(picls, calc_mad, phylota = all_cls)
phylota_table <- data.frame('Cluster.ID' = picls, 'Parent' = prnts,
                            'TaxIDs' = n_taxa, 'GIs' = n_sqs,
                            'Species' = n_species, 'L.min' = mnsqlngs,
                            'L.max' = mxsqlngs, 'MAD' = mads,
                            'Defline.of.longest.sequence' = deflns)
ordr <- order(phylota_table[['TaxIDs']], decreasing = TRUE)
phylota_table <- phylota_table[ordr, ]
write.csv(phylota_table, "../results/Lycium_phylota.csv", row.names = FALSE)

# CLTYPE STATS
# drop clusters of 10
n_taxa <- get_clstr_slot(all_cls, cid=all_cls@cids, slt_nm='ntx')
keep <- names(n_taxa)[n_taxa > 10]
all_cls <- drop_clstrs(all_cls, keep)
table(sapply(all_cls@clstrs@clstrs, function(x) x@typ))

# PLOT
# species treemap
Species_txids <- unique(get_txids(phylota=all_cls,
                                txids=all_cls@txids,
                                rnk='species'))
#remove epithet
Species_txids <- Species_txids[Species_txids != '']

txnms <- get_tx_slot(phylota=all_cls, txid=Species_txids,
                     slt_nm='scnm')
txnms <- sort(txnms, decreasing=TRUE)
p <- plot_phylota_treemap(phylota=all_cls, txids=Species_txids,
                          txnms=txnms, area='nsq', fill='ncl')
png("../results/Lycium_tx_treemap.png", width=2000, height=2000)
print(p + theme(legend.position='none'))
dev.off()
saveRDS(p, "../results/Lycium_tx_treemap.RData")
# cluster treemap
p <- plot_phylota_treemap(phylota=all_cls, cids=all_cls@cids,
                          area='nsq', fill='ntx')
png("../results/Lycium_cluster_treemap.png", width=2000, height=2000)
print(p + theme(legend.position='none'))
dev.off()
saveRDS(p, "../results/Lycium_cluster_treemap.RData")

# REDUCE TO TRIBE
species_only <- drop_by_rank(all_cls, rnk='species', n=10,
                           choose_by= c("pambgs", "age",
                                        "nncltds"),
                           greatest = c(FALSE, FALSE,
                                        TRUE))

# SUMMARISE AND FILTER
# count n species per cid
n_species <- sapply(species_only@cids, function(x) {
  length(unique(get_txids(phylota=species_only, cid=x, rnk='species')))
})
smmry <- summary(species_only)
smmry[['N_taxa']] <- n_species
smmry <- smmry[smmry[['MAD']] > 0.75, ]
smmry <- smmry[order(smmry$N_taxa, decreasing=TRUE), ]

# SELECT
slctd_smmry <- smmry[1:10, ]
slctd_smmry$ID <- as.numeric(slctd_smmry$ID)
slctd <- drop_clstrs(species_only, as.character(slctd_smmry$ID))
write.csv(slctd_smmry, "../results/best_clusters_Lycium.csv")

# # PLOT
# tribes <- unique(get_txids(phylota=all_cls,
#                            txids=all_cls@txids, rnk='tribe'))
# tribes <- tribes[tribes != '']
# txnms <- get_tx_slot(phylota=slctd, txid=tribes, slt_nm='scnm')
# ordd <- order(txnms, decreasing=TRUE)
# p <- plot_phylota_pa(phylota=slctd, cids=slctd@cids,
#                      txids=tribes[ordd], txnms=txnms[ordd])
# png(file.path('results', 'palms_cl_pamap.png'), width=2000, height=2000)
# print(p)
# dev.off()
# saveRDS(p, file.path('figures', 'palms_cl_pamap.RData'))

# OUTPUT
# write out top 10 clusters with most taxa

for(i in seq_along(slctd@cids)) {
  cid <- slctd@cids[i]
  sids <- slctd@clstrs[[cid]]@sids
  txids <- get_txids(slctd, cid=cid, rnk='species')
  scnms <- get_tx_slot(slctd, txids, 'scnm')
  n <- sapply(seq_along(scnms), function(x) 
    sum(scnms[x] == scnms[x:length(scnms)]))
  sq_nm <- paste0(scnms, '_', n)
  infile <- file.path(wd, paste0('sequences', i, '.fasta'))
  write_sqs(phylota=slctd, outfile=infile, sid=sids,
            sq_nm=sq_nm)
}

# ALIGN
#need mafft installed
#for(i in seq_along(slctd@cids)) {
for(i in 1:10) {
  inpt <- file.path(wd, paste0('sequences', i,
                               '.fasta'))
  otpt <- file.path(wd, paste0('alignment', i,'.fasta'))
  system(paste0('mafft --auto ', inpt, ' > ', otpt))
  #muscle
  #system(paste0('muscle -in', inpt, ' -out ', otpt))
}

