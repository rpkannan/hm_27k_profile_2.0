#################### BRCA ########################
setwd("~/Documents/hm_27k_profile_2.0/Output/brca")
brca_n_r = read.csv("normal_rnd_betas.txt",header=TRUE,sep="\t")
barcodes_brca = brca_n_r$Samples
brca_n_r$Samples = NULL
rownames(brca_n_r) = c("B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11","B12","B13","B14","B15","B16","B17","B18","B19","B20","B21","B22","B23","B24","B25","B26","B27")

brca_nr_dist = dist(as.matrix(brca_n_r))

brca_nr_clust = hclust(brca_nr_dist,method="single")
plot(brca_nr_clust,main="Single Linkage")
brca_nr_clust = hclust(brca_nr_dist,method="complete")
plot(brca_nr_clust,main="Complete Linkage")
brca_nr_clust = hclust(brca_nr_dist,method="average")
plot(brca_nr_clust,main="Average Linkage")
brca_nr_clust = hclust(brca_nr_dist,method="ward")
plot(brca_nr_clust,main="Ward Linkage")

###################### COAD ########################
setwd("~/Documents/hm_27k_profile_2.0/Output/coad")
coad_n_r = read.csv("normal_rnd_betas.txt",header=TRUE,sep="\t")
barcodes_coad = coad_n_r$Samples
coad_n_r$Samples = NULL
rownames(coad_n_r) = c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24","C25","C26","C27","C28","C29","C30","C31","C32","C33","C34","C35","C36","C37")

brca_coad = rbind(brca_n_r,coad_n_r)
brca_coad.dist = dist(as.matrix(brca_coad))
brca_coad.dist = dist(as.matrix(brca_coad),method=)

brca_coad_clust = hclust(brca_coad.dist,method="ward")
plot(brca_coad_clust,main="Ward")
brca_coad_clust = hclust(brca_coad.dist,method="single")
plot(brca_coad_clust,main="Single")
brca_coad_clust = hclust(brca_coad.dist,method="complete")
plot(brca_coad_clust,main="complete")