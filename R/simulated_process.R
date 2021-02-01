#######################simulated data
# In order to benchmark the ability of our approach to search and discern the ture DMR or IMR, we resort to a simulation study.
# the sumulated data was constructed from a real WGBS dataset for better capture 
# mouse_new <- read_delim("D:/R/R_project/mouse_new.txt",
#                         "\t", escape_double = FALSE, col_names = FALSE,
#                         trim_ws = TRUE)
# mouse1 <- GenomicRanges::GRanges(seqnames = mouse_new$X1,ranges = IRanges::IRanges(start = mouse_new$X2,end = mouse_new$X3))
# 
# mouseliftOver <- read_delim("D:/R/R_project/mouseliftOver.txt",
#                             "\t", escape_double = FALSE, col_names = FALSE,
#                             trim_ws = TRUE)
# mouse2 <- GenomicRanges::GRanges(seqnames = mouseliftOver$X1,ranges = IRanges::IRanges(start = mouseliftOver$X2,end = mouseliftOver$X3))
# 
# overlaps <- GenomicRanges::findOverlaps(query = mouse2, subject = mouse1,ignore.strand = TRUE)
# length(which(abs((width(mouse2[query]) - width(mouse1[subject])))<300))
# index = which(abs((width(mouse2[query]) - width(mouse1[subject])))<300)
# mouse2[query[index]]
# 取前60个region，30个模拟DMR，30个模拟IMR
# 
# #############simulated data
# example_human_chr22_CpG_report <- read_delim("inst/extdata/example_human_chr22_CpG_report.txt",
#                                 "\t", escape_double = FALSE, col_names = FALSE,
#                                   trim_ws = TRUE)
# 
# human_met <- GenomicRanges::GRanges(seqnames = paste0("chr",example_human_chr22_CpG_report$X1),ranges = IRanges::IRanges(start = example_human_chr22_CpG_report$X2, end = example_human_chr22_CpG_report$X2))
# 
# library(readr)
# example_mouse_chr15_CpG_report <- read_delim("inst/extdata/example_mouse_chr15_CpG_report.txt",
#                             "\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
# 
# mouse_met <- GenomicRanges::GRanges(seqnames = paste0("chr",example_mouse_chr15_CpG_report$X1),ranges = IRanges::IRanges(start = example_mouse_chr15_CpG_report$X2, end = example_mouse_chr15_CpG_report$X2))
# library(readr)
# selected_human_60_CpG_islands <- read_delim("inst/extdata/selected_human_60_CpG_islands.txt",
#             "\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
# 
# library(readr)
# selected_mouse_60_CpG_islands <- read_delim("inst/extdata/selected_mouse_60_CpG_islands.txt",
#             "\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
# 
# human_region <- GenomicRanges::GRanges(seqnames = selected_human_60_CpG_islands$X1,ranges = IRanges::IRanges(start = selected_human_60_CpG_islands$X2,end = selected_human_60_CpG_islands$X3))
# human_region
# mouse_region <- GenomicRanges::GRanges(seqnames = selected_mouse_60_CpG_islands$X1,ranges = IRanges::IRanges(start = selected_mouse_60_CpG_islands$X2,end = selected_mouse_60_CpG_islands$X3))
# human_overlaps <- GenomicRanges::findOverlaps(query = human_met, subject = human_region, ignore.strand = T)
# mouse_overlaps <- GenomicRanges::findOverlaps(query = mouse_met, subject = mouse_region, ignore.strand = T)
# h_query_hits <- S4Vectors::queryHits(human_overlaps)
# h_subj_hits  <- S4Vectors::subjectHits(human_overlaps)
# table(h_subj_hits)
# m_query_hits <- S4Vectors::queryHits(mouse_overlaps)
# m_subj_hits  <- S4Vectors::subjectHits(mouse_overlaps)
# table(m_subj_hits)

# # #simulated 1
# human <- data.frame(seqnames=as.vector(seqnames(human_met[h_query_hits])),site=start(human_met[h_query_hits]),methylated = round(rnorm(NROW(human_met[h_query_hits]),mean = 30,sd=1)), unmethylated = round(rnorm(NROW(human_met[h_query_hits]),mean = 2,sd=2)))
# mouse <- data.frame(seqnames=as.vector(seqnames(mouse_met[m_query_hits])),site=start(mouse_met[m_query_hits]),methylated = round(rnorm(NROW(mouse_met[m_query_hits]),mean = 30,sd=1)), unmethylated = round(rnorm(NROW(mouse_met[m_query_hits]),mean = 2,sd=2)))
# names(human)[names(human) == "site"] <- 'X2'
# names(mouse)[names(mouse) == "site"] <- 'X2'
# human_merage <- merge(human,example_human_chr22_CpG_report,all=TRUE)
# mouse_merage <- merge(mouse,example_mouse_chr15_CpG_report,all=TRUE)
# mouse_merage[which(mouse_merage$unmethylated<0),]$unmethylated = 0
# human_merage[which(human_merage$unmethylated<0),]$unmethylated = 0
# write.table(human_merage,file = '../human_merge.txt',row.names = F,col.names = F,sep = "\t")
# write.table(mouse_merage,file = '../mouse_merge.txt',row.names = F,col.names = F,sep = "\t")
# # sed -i 's/"//g' human_merge.txt
# # sed -i 's/"//g' mouse_merge.txt
# # awk -F "\t" '{if($3=="NA"){print $5FS$1FS$6FS$7FS$8FS$9FS$10} else {print $5FS$1FS$6FS$3FS$4FS$9FS$10}}' human_merge.txt > human_simulated_1.txt &
# # awk -F "\t" '{if($3=="NA"){print $5FS$1FS$6FS$7FS$8FS$9FS$10} else {print $5FS$1FS$6FS$3FS$4FS$9FS$10}}' mouse_merge.txt > mouse_simulated_1.txt &

#
# #simulated 2
# human <- data.frame(seqnames=as.vector(seqnames(human_met[h_query_hits])),site=start(human_met[h_query_hits]),methylated = round(rnorm(NROW(human_met[h_query_hits]),mean = 30, sd=1)), unmethylated = round(rnorm(NROW(human_met[h_query_hits]),mean = 2,sd=2)))
# mouse <- data.frame(seqnames=as.vector(seqnames(mouse_met[m_query_hits])),site=start(mouse_met[m_query_hits]),methylated = round(rnorm(NROW(mouse_met[m_query_hits]),mean = 2, sd=1)), unmethylated = round(rnorm(NROW(mouse_met[m_query_hits]),mean = 30,sd=2)))
# names(human)[names(human) == "site"] <- 'X2'
# names(mouse)[names(mouse) == "site"] <- 'X2'
# human_merage <- merge(human,example_human_chr22_CpG_report,all=TRUE)
# mouse_merage <- merge(mouse,example_mouse_chr15_CpG_report,all=TRUE)
# mouse_merage[which(mouse_merage$methylated<0),]$methylated = 0
# human_merage[which(human_merage$unmethylated<0),]$unmethylated = 0
# write.table(human_merage,file = '../human_merge.txt',row.names = F,col.names = F,sep = "\t")
# write.table(mouse_merage,file = '../mouse_merge.txt',row.names = F,col.names = F,sep = "\t")
# # sed -i 's/"//g' human_merge.txt
# # sed -i 's/"//g' mouse_merge.txt
# # awk -F "\t" '{if($3=="NA"){print $5FS$1FS$6FS$7FS$8FS$9FS$10} else {print $5FS$1FS$6FS$3FS$4FS$9FS$10}}' human_merge.txt > human_simulated_2.txt &
# # awk -F "\t" '{if($3=="NA"){print $5FS$1FS$6FS$7FS$8FS$9FS$10} else {print $5FS$1FS$6FS$3FS$4FS$9FS$10}}' mouse_merge.txt > mouse_simulated_2.txt &




# #simulated 3
# human <- data.frame(seqnames = character(),site= integer(),methylated=numeric(),unmethylated=numeric())
# for (i in unique(h_subj_hits)) {
#   a = h_query_hits[which(h_subj_hits==i)]
#   b = data.frame(seqnames = as.vector(seqnames(human_met[a])),site=start(human_met[a]),methylated=rep(1:length(a)),unmethylated=rep(length(a):1))
#   human <- rbind(human,b)
# }
# mouse <- data.frame(seqnames = character(),site= integer(),methylated=numeric(),unmethylated=numeric())
# for (i in unique(m_subj_hits)) {
#   a = m_query_hits[which(m_subj_hits==i)]
#   b = data.frame(seqnames = as.vector(seqnames(mouse_met[a])),site=start(mouse_met[a]),methylated=rep(length(a):1),unmethylated=rep(1:length(a)))
#   mouse <- rbind(mouse,b)
# }
# names(human)[names(human) == "site"] <- 'X2'
# names(mouse)[names(mouse) == "site"] <- 'X2'
# human_merage <- merge(human,example_human_chr22_CpG_report,all=TRUE)
# mouse_merage <- merge(mouse,example_mouse_chr15_CpG_report,all=TRUE)
# write.table(mouse_merage,file = '../mouse_merge.txt',row.names = F,col.names = F,sep = "\t")
# # sed -i 's/"//g' human_merge.txt
# # sed -i 's/"//g' mouse_merge.txt
# # awk -F "\t" '{if($3=="NA"){print $5FS$1FS$6FS$7FS$8FS$9FS$10} else {print $5FS$1FS$6FS$3FS$4FS$9FS$10}}' human_merge.txt > human_simulated_3.txt &
# # awk -F "\t" '{if($3=="NA"){print $5FS$1FS$6FS$7FS$8FS$9FS$10} else {print $5FS$1FS$6FS$3FS$4FS$9FS$10}}' mouse_merge.txt > mouse_simulated_3.txt &
#
# #
# #step1. read files
# file <- system.file("extdata", "human_simulated_2.txt", package = "BSDMR")
# human_met <- read_methylation_report(file,min_coverage=10,type = "CG")
# file <- system.file("extdata", "mouse_simulated_2.txt", package = "BSDMR")
# mouse_met <- read_methylation_report(file,min_coverage=10,type = "CG")
# ##read annotation file (customized region)
# file <- system.file("extdata", "human_chr22_annotation.txt", package = "BSDMR")
# human_anno <- read_annotation(file)
# file <- system.file("extdata", "mouse_chr15_annotation.txt", package = "BSDMR")
# mouse_anno <- read_annotation(file)
# 
# #step2. cluster_sites_to_region
# # human_region <- cluster_sites_to_region(methylation = human_met,
# #                                         annotation = human_anno,
# #                                         min_sites_number = 10,
# #                                         max_distance_between_sites = 100,
# #                                         is_parallel = TRUE)
# #
# # #step3. liftOver map the human_region to the annotated region of another species
# # filePath <- system.file("extdata", "hg38ToMm10.over.chain", package = "BSDMR")
# # mouse_region <- change_genomic_coordinate(human_region,filePath,mouse_anno)
# 
# # # or verse
# mouse_region <- cluster_sites_to_region(methylation = mouse_met,
#                                         annotation = mouse_anno,
#                                         min_sites_number = 10,
#                                         max_distance_between_sites = 100,
#                                         is_parallel = TRUE)
# filePath <- system.file("extdata", "mm10ToHg38.over.chain", package = "BSDMR")
# human_region <- change_genomic_coordinate(mouse_region,filePath,human_anno)
# 
# #step4. create region object as the input of modeling
# human_obj <- create_region_object(human_met, human_region)
# mouse_obj <- create_region_object(mouse_met, mouse_region)
# 
# #step5. infer profiles of reion object
# human_basis_profile <- create_rbf_object(M = 8)
# human_fit_profiles <- infer_profiles_vb(X = human_obj$met, basis = human_basis_profile,
#                                         is_parallel = TRUE, vb_max_iter = 100)
# 
# human_basis_mean <- create_rbf_object(M = 0)
# human_fit_mean <- infer_profiles_vb(X = human_obj$met, basis = human_basis_mean,
#                                     is_parallel = TRUE, vb_max_iter = 100)
# 
# mouse_basis_profile <- create_rbf_object(M = 8)
# mouse_fit_profiles <- infer_profiles_vb(X = mouse_obj$met, basis = mouse_basis_profile,
#                                         is_parallel = TRUE, vb_max_iter = 100)
# 
# mouse_basis_mean <- create_rbf_object(M = 0)
# mouse_fit_mean <- infer_profiles_vb(X = mouse_obj$met, basis = mouse_basis_mean,
#                                     is_parallel = TRUE, vb_max_iter = 100)
# 
# #step6. computer the adjusted consine distance of profiles and measure the similarity
# similarity <- adjusted_cosine_similarity(queryProfiles=human_fit_profiles,subjectProfiles=mouse_fit_profiles)
# which(similarity>0.9)
# length(similarity[which(!is.na(similarity))])
# i=23
# plot_infer_profiles(region = i, obj_prof = human_fit_profiles,obj_mean = human_fit_mean,
#                     obs = human_obj, title = paste0("Gene ID ",human_obj$anno$id[i]))
# plot_infer_profiles(region = i, obj_prof = mouse_fit_profiles,obj_mean = mouse_fit_mean,
#                     obs = mouse_obj, title = paste0("Gene ID ",mouse_obj$anno$id[i]))
# 
# library(readr)
# selected_human_60_CpG_islands <- read_delim("inst/extdata/selected_human_60_CpG_islands.txt",
#             "\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
# 
# library(readr)
# selected_mouse_60_CpG_islands <- read_delim("inst/extdata/selected_mouse_60_CpG_islands.txt",
#             "\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)
# 
# mouse_region_60 <- GenomicRanges::GRanges(seqnames = selected_mouse_60_CpG_islands$X1,ranges = IRanges::IRanges(start = selected_mouse_60_CpG_islands$X2,end = selected_mouse_60_CpG_islands$X3))
# human_region_60 <- GenomicRanges::GRanges(seqnames = selected_human_60_CpG_islands$X1,ranges = IRanges::IRanges(start = selected_human_60_CpG_islands$X2,end = selected_human_60_CpG_islands$X3))
# overlaps1 <- GenomicRanges::findOverlaps(query = human_region, subject = human_region_60, ignore.strand = T)
# length(unique(subjectHits(overlaps1)))
# overlaps2 <- GenomicRanges::findOverlaps(query = mouse_region, subject = mouse_region_60, ignore.strand = T)
# length(unique(subjectHits(overlaps2)))
# #查看mouse_reion中成功liftOver中转换成功的区域有多少
# length(which(mouse_region$id!="liftOver FAILED"))
# #查看有多少DMR被找出来了
# length(which(similarity>0.9))
# #查看最终多少DMR是真的DMR
# overlaps3 <- GenomicRanges::findOverlaps(query = mouse_region[which(similarity>0.9)], subject = mouse_region_60, ignore.strand = T)
# overlaps4 <- GenomicRanges::findOverlaps(query = human_region[which(similarity>0.9)], subject = human_region_60, ignore.strand = T)
# #length(unique(subjectHits(overlaps3)))
# #length(unique(subjectHits(overlaps4)))
# length(subjectHits(overlaps3))
# length(subjectHits(overlaps4))
# 


