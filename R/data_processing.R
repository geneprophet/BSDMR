# mouse_met <- read_methylation_report("/data40T/kanghe/WGBS/mouse_C57BL/mouse_14",min_coverage=4)
# mouse_anno <- read_annotation("/data40T/kanghe/R/Mouse.gtf")
# human_anno <- read_annotation("/data40T/kanghe/R/Human.gtf")
# human_met <- read_methylation_report("/data40T/kanghe/WGBS/human/human_10_CX_report.txt",min_coverage=4)
# human_region <- cluster_sites_to_region(methylation = human_met,annotation = human_anno,is_parallel = TRUE)
# filePath <- system.file("extdata", "hg38ToMm10.over.chain", package = "BSDMR")
# mouse_region <- change_genomic_coordinate(human_region,filePath,mouse_anno)
# human_basis_profile <- create_rbf_object(M = 8)
# human_basis_mean <- create_rbf_object(M = 0)
# human_obj <- create_region_object(human_met, human_region)
# human_fit_profiles <- infer_profiles_vb(X = human_obj$met, model = "binomial",basis = human_basis_profile, is_parallel = TRUE, vb_max_iter = 100)
# human_fit_mean <- infer_profiles_vb(X = human_obj$met, model = "binomial",basis = human_basis_mean, is_parallel = TRUE, vb_max_iter = 100)
# 
# mouse_basis_profile <- create_rbf_object(M = 8)
# mouse_basis_mean <- create_rbf_object(M = 0)
# mouse_obj <- create_region_object(mouse_met, mouse_region)
# mouse_fit_profiles <- infer_profiles_vb(X = mouse_obj$met, model = "binomial",basis = mouse_basis_profile, is_parallel = TRUE, vb_max_iter = 100)
# mouse_fit_mean <- infer_profiles_vb(X = mouse_obj$met, model = "binomial",basis = mouse_basis_mean, is_parallel = TRUE, vb_max_iter = 100)




