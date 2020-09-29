## code to prepare `DATASET` dataset goes here

# usethis::use_data(DATASET, overwrite = TRUE)

generate_human_met <- function(){
  file <- system.file("extdata", "example_human_CpG_report.txt", package = "BSDMR")
  human_met <- read_methylation_report(file,min_coverage=4)
  return(human_met)
}
generate_mouse_met <- function(){
  file <- system.file("extdata", "example_mouse_CpG_report.txt", package = "BSDMR")
  mouse_met <- read_methylation_report(file,min_coverage=4)
  return(mouse_met)
}

generate_human_anno <- function(){
  file <- system.file("extdata", "example_human_annotation.txt", package = "BSDMR")
  human_anno <- read_annotation(file)
  return(human_anno)
}

generate_mouse_anno <- function(){
  file <- system.file("extdata", "example_mouse_annotation.txt", package = "BSDMR")
  mouse_anno <- read_annotation(file)
  return(mouse_anno)
}

generate_human_region <- function(){
  human_region <- cluster_sites_to_region(methylation = human_met,annotation = human_anno,is_parallel = TRUE)
  return(human_region)
}
generate_mouse_region <- function(){
   filePath <- system.file("extdata", "hg38ToMm10.over.chain", package = "BSDMR")
   mouse_region <- change_genomic_coordinate(human_region,filePath,mouse_anno)
}
generate_human_obj <- function(){
  human_obj <- create_region_object(human_met, human_region)
}
generate_mouse_obj <- function(){
  mouse_obj <- create_region_object(mouse_met, mouse_region)
}
generate_human_fit_profiles <- function(){
   human_basis_profile <- create_rbf_object(M = 8)
   human_obj <- create_region_object(human_met, human_region)
   human_fit_profiles <- infer_profiles_vb(X = human_obj$met, model = "binomial",
      basis = human_basis_profile, is_parallel = TRUE, vb_max_iter = 100)
}
generate_human_fit_mean <- function(){
  human_basis_mean <- create_rbf_object(M = 0)
  human_obj <- create_region_object(human_met, human_region)
  human_fit_mean <- infer_profiles_vb(X = human_obj$met, model = "binomial",
                                      basis = human_basis_mean, is_parallel = TRUE, vb_max_iter = 100)
}
generate_mouse_fit_profiles <- function(){
  mouse_basis_profile <- create_rbf_object(M = 8)
  mouse_obj <- create_region_object(mouse_met, mouse_region)
  mouse_fit_profiles <- infer_profiles_vb(X = mouse_obj$met, model = "binomial",
                                          basis = mouse_basis_profile, is_parallel = TRUE, vb_max_iter = 100)
}
generate_mouse_fit_mean <- function(){
  mouse_basis_mean <- create_rbf_object(M = 0)
  mouse_obj <- create_region_object(mouse_met, mouse_region)
  mouse_fit_mean <- infer_profiles_vb(X = mouse_obj$met, model = "binomial",
                                      basis = mouse_basis_mean, is_parallel = TRUE, vb_max_iter = 100)
}

generate_similarity <- function(){
  similarity <- adjusted_cosine_similarity(queryProfiles=human_fit_profiles,subjectProfiles=mouse_fit_profiles)
}
