##3k heart human and mouse
# mouse_anno <- read_annotation("/data40T/kanghe/R/MOUSE_3K.gtf")
# human_anno <- read_annotation("/data40T/kanghe/R/HUMAN_3K.gtf")
# human_met_heart <- read_methylation_report("/data40T/kanghe/WGBS/human/heart/heart.txt",min_coverage=4)
# mouse_met_heart <- read_methylation_report("/data40T/kanghe/WGBS/mouse_C57BL/heart/heart.txt",min_coverage=4)
# human_region_heart <- cluster_sites_to_region(methylation = human_met_heart, annotation = human_anno, is_parallel = TRUE)
# filePath <- system.file("extdata", "hg38ToMm10.over.chain", package = "BSDMR")
# mouse_region_heart <- change_genomic_coordinate(human_region_heart,filePath,mouse_anno)
# human_obj_heart <- create_region_object(human_met_heart, human_region_heart)
# mouse_obj_heart <- create_region_object(mouse_met_heart, mouse_region_heart)
# human_basis_profile <- create_rbf_object(M = 8)
# human_basis_mean <- create_rbf_object(M = 0)
# human_fit_profiles_heart <- infer_profiles_vb(X = human_obj_heart$met, model = "binomial",basis = human_basis_profile, is_parallel = TRUE, vb_max_iter = 100)
# human_fit_mean_heart <- infer_profiles_vb(X = human_obj_heart$met, model = "binomial",basis = human_basis_mean, is_parallel = TRUE, vb_max_iter = 100)
# 
# mouse_basis_profile <- create_rbf_object(M = 8)
# mouse_basis_mean <- create_rbf_object(M = 0)
# mouse_fit_profiles_heart <- infer_profiles_vb(X = mouse_obj_heart$met, model = "binomial",basis = mouse_basis_profile, is_parallel = TRUE, vb_max_iter = 100)
# mouse_fit_mean_heart <- infer_profiles_vb(X = mouse_obj_heart$met, model = "binomial",basis = mouse_basis_mean, is_parallel = TRUE, vb_max_iter = 100)
# similarity_heart <- adjusted_cosine_similarity(queryProfiles=human_fit_profiles_heart,subjectProfiles=mouse_fit_profiles_heart)

# #查看心脏组织的情况
# load('../heart_30.Rdata')
# h_dmr = human_region_heart[which(similarity_heart<0.2)]
# m_dmr = mouse_region_heart[which(similarity_heart<0.2)]
# 
# h_imr = human_region_heart[which(similarity_heart>0.9)]
# m_imr = mouse_region_heart[which(similarity_heart>0.9)]
# 
# covplot(h_dmr)
# txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
# h_dmr_anno = annotatePeak(h_dmr,tssRegion = c(-3000,3000),TxDb = txdb,annoDb = "org.Hs.eg.db")
# plotAnnoPie(h_dmr_anno)
# vennpie(h_dmr_anno)
# upsetplot(h_dmr_anno)
# 
# covplot(h_imr)
# txdb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
# h_imr_anno = annotatePeak(h_imr,tssRegion = c(-3000,3000),TxDb = txdb,annoDb = "org.Hs.eg.db")
# plotAnnoPie(h_imr_anno)
# vennpie(h_imr_anno)
# upsetplot(h_imr_anno)