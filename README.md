
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BSDMR

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/geneprophet/BSDMR?branch=master&svg=true)](https://ci.appveyor.com/project/geneprophet/BSDMR)
<!-- badges: end -->

The goal of BSDMR is to find DMRs or CMRs by high order methylation
profiles between species.

## Installation

You can install the development version from
[GitHub](https://github.com/geneprophet/BSDMR) with:

``` r
# install.packages("devtools")
devtools::install_github("geneprophet/BSDMR")
```

## Example

This is a basic example which shows you how to use the package:

``` r
#library the package
library(BSDMR)

#example pipeline

#step1. read files
##read methylation report (the result of bismark)
file <- system.file("extdata", "example_human_chr22_CpG_report.txt", package = "BSDMR")
human_met <- read_methylation_report(file,min_coverage=4,type = "CG")
#> reading methylation report file ...
file <- system.file("extdata", "example_mouse_chr15_CpG_report.txt", package = "BSDMR")
mouse_met <- read_methylation_report(file,min_coverage=4,type = "CG")
#> reading methylation report file ...

##read annotation file (customized region)
file <- system.file("extdata", "human_chr22_annotation.txt", package = "BSDMR")
human_anno <- read_annotation(file)
#> reading annotation file...
file <- system.file("extdata", "mouse_chr15_annotation.txt", package = "BSDMR")
mouse_anno <- read_annotation(file)
#> reading annotation file...

#step2. cluster_sites_to_region
human_region <- cluster_sites_to_region(methylation = human_met,
                                        annotation = human_anno,
                                        min_sites_number = 10,
                                        max_distance_between_sites = 200,
                                        is_parallel = TRUE)

#step3. liftOver map the human_region to the annotated region of another species(mouse)
filePath <- system.file("extdata", "hg38ToMm10.over.chain", package = "BSDMR")
mouse_region <- change_genomic_coordinate(human_region,filePath,mouse_anno)

#step4. create region object as the input of modeling
human_obj <- create_region_object(human_met, human_region)
#> Creating methylation regions ...
mouse_obj <- create_region_object(mouse_met, mouse_region)
#> Creating methylation regions ...

#step5. infer profiles of reion object
human_basis_profile <- create_rbf_object(M = 8)
human_fit_profiles <- infer_profiles_vb(X = human_obj$met, basis = human_basis_profile, 
                                        is_parallel = TRUE, vb_max_iter = 100)
#> infering the methylation profiles by variational bayes....

human_basis_mean <- create_rbf_object(M = 0)
human_fit_mean <- infer_profiles_vb(X = human_obj$met, basis = human_basis_mean,
                                    is_parallel = TRUE, vb_max_iter = 100)
#> infering the methylation profiles by variational bayes....

mouse_basis_profile <- create_rbf_object(M = 8)
mouse_fit_profiles <- infer_profiles_vb(X = mouse_obj$met, basis = mouse_basis_profile, 
                                        is_parallel = TRUE, vb_max_iter = 100)
#> infering the methylation profiles by variational bayes....

mouse_basis_mean <- create_rbf_object(M = 0)
mouse_fit_mean <- infer_profiles_vb(X = mouse_obj$met, basis = mouse_basis_mean, 
                                    is_parallel = TRUE, vb_max_iter = 100)
#> infering the methylation profiles by variational bayes....

#step6. computer the adjusted consine distance of profiles and measure the similarity
similarity <- adjusted_cosine_similarity(queryProfiles = human_fit_profiles,subjectProfiles = mouse_fit_profiles)

#step7. visualization
which(similarity>0.9)
#>  [1] 160 170 419 439 442 483 525 535 565 573 577 644 651 679 744 771 801 831 838
#> [20] 873 971 989
i=160
plot_infer_profiles(region = i, obj_prof = human_fit_profiles,obj_mean = human_fit_mean,
                    obs = human_obj, title = paste0("Gene ID ",human_obj$anno$id[i]))
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
plot_infer_profiles(region = i, obj_prof = mouse_fit_profiles,obj_mean = mouse_fit_mean,
                    obs = mouse_obj, title = paste0("Gene ID ",mouse_obj$anno$id[i]))
```

<img src="man/figures/README-example-2.png" width="100%" />
