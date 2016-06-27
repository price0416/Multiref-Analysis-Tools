
#A few quick commands to merge our DE output files with our gene map and remove any duplicates.

#Change the header of the original reference gene to 'Gene' before proceeding.  
mapfile <- read.table("CMCP6_YJ016_gene_mapping.dat", sep="\t", header=TRUE)
de1_file <- read.table(file="YJ016_HS_vs_CMCP6ref_ASW_difExp_results.csv", header=TRUE, sep=",")
table1 <- merge(mapfile, de1_file, by="Gene")
tmpTable <- table1
table1<-subset(tmpTable, !duplicated(Gene))
write.csv(table1, file="YJ016_HS_vs_CMCP6ref_ASW_difExp_results_mapMerged.csv")