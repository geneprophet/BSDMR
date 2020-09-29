#' @keywords internal
#' @title BSDMR : Find DMRs between species by high order methylation profile
#' 
#' @description Extracting higher order methylation features by GLM to find DMRs between species.
#'              A region-based method to find DMR.
#' 
#' @name BSDMR
#' 
#' @author Hongen Kang  \email{geneprophet@163.com}
#' 
#' 
#' @rawNamespace importFrom(magrittr,"%>%")
#' @rawNamespace importFrom(data.table,":=")
#' @import GenomicRanges ggplot2 foreach
#' @importFrom assertthat assert_that
#' @importFrom stats pnorm dbinom dnorm sd
#' @importFrom Rcpp evalCpp
#' @importFrom readr read_delim cols col_character col_integer
#' @importFrom S4Vectors queryHits 
#' @importFrom IRanges IRanges
#' @importFrom rtracklayer liftOver import.chain
#' @importFrom BiocGenerics as.vector unique strand unique
#' @importFrom utils head tail
#' @importFrom stats cor
#' @importFrom matrixcalc matrix.trace
#' 
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @useDynLib BSDMR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

if (nzchar(chk) && chk == "TRUE") {
  # use 2 cores in CRAN/Travis/AppVeyor
  num_workers <- 2L
} else {
  # use all cores in devtools::test()
  num_workers <- parallel::detectCores()
}
.datatable.aware <- TRUE
.onLoad <- function(libname = find.package("BSDMR"), pkgname = "BSDMR"){
  # CRAN Note avoidance
  if (getRversion() >= "2.15.1")
    utils::globalVariables(
      # sample file names from taxstats
      c(# we use the magrittr pipe
        ".",
        "i",
        "X4",
        "X5",
        "X6"
      )
    )
  invisible()
}
# 
# #############################################################
# #example pipeline
# #step1. read files
# ##read methylation report (the result of bismark)
# file <- system.file("extdata", "example_human_CpG_report.txt", package = "BSDMR")
# human_met <- read_methylation_report(file,min_coverage=4,type = "CG")
# file <- system.file("extdata", "example_mouse_CpG_report.txt", package = "BSDMR")
# mouse_met <- read_methylation_report(file,min_coverage=4,type = "CG")
# file <- system.file("extdata", "example_human_annotation.txt", package = "BSDMR")
# human_anno <- read_annotation(file)
# ##read annotation file (customized region)
# file <- system.file("extdata", "example_human_annotation.txt", package = "BSDMR")
# human_anno <- read_annotation(file)
# file <- system.file("extdata", "example_mouse_annotation.txt", package = "BSDMR")
# mouse_anno <- read_annotation(file)
# 
# #step2. cluster_sites_to_region
# human_region <- cluster_sites_to_region(methylation = human_met,
#                                         annotation = human_anno,
#                                         min_sites_number = 10,
#                                         max_distance_between_sites = 200,
#                                         is_parallel = TRUE)
# 
# #step3. liftOver map the human_region to the annotated region of another species
# filePath <- system.file("extdata", "hg38ToMm10.over.chain", package = "BSDMR")
# targetGRanges <- change_genomic_coordinate(human_region,filePath,mouse_anno)
# 
# #step4. create region object as the input of modeling
# human_obj <- create_region_object(human_met, human_anno)
# mouse_obj <- create_region_object(mouse_met, targetGRanges)
# 
# #step5. infer profiles of reion object
# human_basis_profile <- create_rbf_object(M = 8)
# human_fit_profiles <- infer_profiles_vb(X = human_obj$met, model = "binomial",
#                                         basis = human_basis_profile, is_parallel = TRUE, vb_max_iter = 100)
# 
# human_basis_mean <- create_rbf_object(M = 0)
# human_fit_mean <- infer_profiles_vb(X = human_obj$met, model = "binomial",
#                                     basis = human_basis_mean, is_parallel = TRUE, vb_max_iter = 100)
# 
# mouse_basis_profile <- create_rbf_object(M = 8)
# mouse_fit_profiles <- infer_profiles_vb(X = mouse_obj$met, model = "binomial",
#                                         basis = mouse_basis_profile, is_parallel = TRUE, vb_max_iter = 100)
# 
# mouse_basis_mean <- create_rbf_object(M = 0)
# mouse_fit_mean <- infer_profiles_vb(X = mouse_obj$met, model = "binomial",
#                                     basis = mouse_basis_mean, is_parallel = TRUE, vb_max_iter = 100)
# 
# #step6. computer the adjusted consine distance of profiles and measure the similarity
# similarity <- adjusted_cosine_similarity(queryProfiles=human_fit_profiles,queryObj=human_obj,subjectProfiles=mouse_fit_profiles,subjectObj=mouse_obj)
# 
# #step7. visualization
# which(similarity$similarity<0.2)
# similarity[43,]
# plot_infer_profiles(region = 43, obj_prof = human_fit_profiles,obj_mean = human_fit_mean,
#                     obs = human_obj, title = paste0("Gene ID ",human_obj$anno$id[43]))
# plot_infer_profiles(region = 4, obj_prof = mouse_fit_profiles,obj_mean = mouse_fit_mean,
#                     obs = mouse_obj, title = paste0("Gene ID ",mouse_obj$anno$id[4]))
