###3-3. Functional module profiling with the curated gut metabolic module (GMMs) using Omixer-RPM
# This step was performed in the R environment using the R packages, omixer-rpmR, based on KEGG feature tables as input.

library(omixerRpm)
pwy_module <- read.table("/path/to/humann_genefamilies_kegg_merged_relab.txt", header = T, row.names = NULL, sep = "\t")

db_gmm <- loadDB("GMMs.v1.07")
result_gmm <- rpm(pwy_module, module.db = db_gmm, minimum.coverage=0.666, annotation = 1)
result_gmm_abn <- asDataFrame(result_gmm, "abundance")
write.table(result_gmm_abn, "/path/to/GMM_abundance.txt", sep="\t", row.names=FALSE)