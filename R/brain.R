############################################
#3k brain human and mouse
# mouse_anno <- read_annotation("/data40T/kanghe/R/MOUSE_3K.gtf")
# human_anno <- read_annotation("/data40T/kanghe/R/HUMAN_3K.gtf")
# human_met_brain <- read_methylation_report("/data40T/kanghe/WGBS/human/brain/human_12_brain",min_coverage=4,type = "CHG")
# mouse_met_brain <- read_methylation_report("/data40T/kanghe/WGBS/mouse_C57BL/brain/mouse_6_brain",min_coverage=4,type = "CHG")
# human_region_brain <- cluster_sites_to_region(methylation = human_met_brain, annotation = human_anno, is_parallel = TRUE)
# filePath <- system.file("extdata", "hg38ToMm10.over.chain", package = "BSDMR")
# mouse_region_brain <- change_genomic_coordinate(human_region_brain,filePath,mouse_anno)
# human_obj_brain <- create_region_object(human_met_brain, human_region_brain)
# mouse_obj_brain <- create_region_object(mouse_met_brain, mouse_region_brain)
# human_basis_profile <- create_rbf_object(M = 8)
# human_basis_mean <- create_rbf_object(M = 0)
# human_fit_profiles_brain <- infer_profiles_vb(X = human_obj_brain$met, model = "binomial",basis = human_basis_profile, is_parallel = TRUE, vb_max_iter = 100)
# human_fit_mean_brain <- infer_profiles_vb(X = human_obj_brain$met, model = "binomial",basis = human_basis_mean, is_parallel = TRUE, vb_max_iter = 100)
# 
# mouse_basis_profile <- create_rbf_object(M = 8)
# mouse_basis_mean <- create_rbf_object(M = 0)
# mouse_fit_profiles_brain <- infer_profiles_vb(X = mouse_obj_brain$met, model = "binomial",basis = mouse_basis_profile, is_parallel = TRUE, vb_max_iter = 100)
# mouse_fit_mean_brain <- infer_profiles_vb(X = mouse_obj_brain$met, model = "binomial",basis = mouse_basis_mean, is_parallel = TRUE, vb_max_iter = 100)
# similarity_brain <- adjusted_cosine_similarity(queryProfiles=human_fit_profiles_brain,subjectProfiles=mouse_fit_profiles_brain)
# 

# #查看脑组织的情况
# load('../brain_30.Rdata')
# h_dmr = human_region_brain[which(similarity_brain<0.2)]
# m_dmr = mouse_region_brain[which(similarity_brain<0.2)]
# 
# h_imr = human_region_brain[which(similarity_brain>0.9)]
# m_imr = mouse_region_brain[which(similarity_brain>0.9)]
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