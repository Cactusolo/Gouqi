
genSumTable <- function(clstrs_obj) {
  sqlngs <- getClMdLn(clstrs_obj, clstrs_obj@clstr_ids)
  n_taxa <- getClNTx(clstrs_obj, clstrs_obj@clstr_ids)
  sq_dfs <- rep(NA, length(n_taxa))
  for(i in seq_along(clstrs_obj@clstr_ids)) {
    id <- clstrs_obj@clstr_ids[i]
    sq_dfs[i] <- getSqDfs(clstrs_obj, id=id, prse=0.1)
  }
  names(sq_dfs) <- clstrs_obj@clstr_ids
  top_ids <- names(sort(n_taxa, decreasing=TRUE))
  cis <- sapply(top_ids, function(x) clstrs_obj@clstrs[[x]][['ci']])
  ntids <- sapply(top_ids, function(x) clstrs_obj@clstrs[[x]][['n_ti']])
  seeds <- sapply(top_ids, function(x) clstrs_obj@clstrs[[x]][['seed_gi']])
  smmry <- data.frame('id'=cis, 'n_tx'=n_taxa[top_ids],
                      'med_sq_lngth'=sqlngs[top_ids],
                      'mad'=mad_scrs[top_ids], 'seed'=seeds[top_ids],
                      'sq_df'=sq_dfs[top_ids])
  row.names(smmry) <- NULL
  smmry
}

rmFls <- function(fls) {
  if(length(fls) > 0) {
    for(fl in fls) {
      file.remove(fl)
    }
  }
}
