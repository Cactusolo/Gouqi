library("rotl")
library("Taxonstand")

query <- read.csv("Gouqi_update.csv", header=FALSE, stringsAsFactors=FALSE)

#OTT
resolved_names <- tnrs_match_names(query$V1, context_name = "Land plants")
write.csv(resolved_names, "Gouqi_checkname_OttQery.csv", row.names=FALSE)

#TPL
TPL_results <- TPL(query$V1, infra = TRUE, diffchar=2, max.distance=1, corr=TRUE, encoding = "UTF-8")
write.csv(TPL_results, "Gouqi_checkname_TPLQery.csv", row.names=FALSE)
