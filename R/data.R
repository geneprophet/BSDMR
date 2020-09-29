#' @title a example output of \code{read_methylation_report}
#'
#' @description a samll subset methylation data of human chr22
#'
#' @format A \code{GRanges} object.
#'
#'   The GRanges object contains three additional metadata columns:
#'   \itemize{ \item{ \code{methylated_reads}: the number if methylated reads for each C sites.}  
#'   \item{ \code{total_reads}: Total number of reads for each C sites.} 
#'   \item{ \code{type}:the C-context type: CG,CHG, or CHH} 
#'   }
#'   
#' @return  methylation data
#'
#' @seealso \code{\link{read_methylation_report}}
#'   
"human_met"

#' @title a example output of \code{read_methylation_report}
#'
#' @description a samll subset methylation data of mouse chr15
#'
#' @format A \code{GRanges} object.
#'
#'   The GRanges object contains three additional metadata columns:
#'   \itemize{ \item{ \code{methylated_reads}: the number if methylated reads for each C sites.}  
#'   \item{ \code{total_reads}: Total number of reads for each C sites.} 
#'   \item{ \code{type}:the C-context type: CG,CHG, or CHH} 
#'   }
#' 
#' @return  methylation data
#'
#' @seealso \code{\link{read_methylation_report}}
#'   
"mouse_met"

#' @title a example output of \code{read_annotation}
#'
#' @description a samll subset human genomic annotation of chr22
#'
#' @format A \code{GRanges} object.
#'
#' The GRanges object contains three additional metadata columns:
#'      \itemize{ \item{ \code{id}: EnsemblID}  
#'      \item{ \code{center}: the center of the region,(strat+end)/2} 
#'      \item{ \code{annotation}:the annotation of the region, it can be promoter,downstream,gene,UTR...} }
#' 
#' 
#' @return  annotation data
#'
#' @seealso \code{\link{read_annotation}}
#'   
"human_anno"


#' @title a example output of \code{read_annotation}
#'
#' @description a samll subset mouse genomic annotation of chr15
#'
#' @format A \code{GRanges} object.
#'
#' The GRanges object contains three additional metadata columns:
#'      \itemize{ \item{ \code{id}: EnsemblID}  
#'      \item{ \code{center}: the center of the region,(strat+end)/2} 
#'      \item{ \code{annotation}:the annotation of the region, it can be promoter,downstream,gene,UTR...} }
#' 
#' 
#' @return  annotation data
#'
#' @seealso \code{\link{read_annotation}}
#'   
"mouse_anno"


#' @title a example output of \code{cluster_sites_to_region }
#'
#' @description \code{human_region <- cluster_sites_to_region(methylation = human_met,
#'                                annotation = human_anno,is_parallel = TRUE)}
#'
#' @format  A \code{GRanges} object.
#' The GRanges object contains three additional metadata columns:
#'      \itemize{ \item{ \code{id}: EnsemblID}  
#'      \item{ \code{center}: the center of the region,(strat+end)/2} 
#'      \item{ \code{annotation}:the annotation of the region, it can be promoter,downstream,gene,UTR...} }
#'
#'     
#' @return annotation data
#'
#' @seealso \code{\link{cluster_sites_to_region}}
#' 
"human_region"



#' @title a example output of \code{change_genomic_coordinate }
#'
#' @description \code{ filePath <- system.file("extdata", "hg38ToMm10.over.chain", package = "BSDMR")
#'                    mouse_region <- change_genomic_coordinate(human_region,filePath,mouse_anno)} 
#'
#' @format  A \code{GRanges} object.
#' The GRanges object contains three additional metadata columns:
#'      \itemize{ \item{ \code{id}: EnsemblID or liftOver FAILED,indicates the region mapping failed between seocies}  
#'      \item{ \code{center}: the center of the region,(strat+end)/2, or NA if liftOver FAILED} 
#'      \item{ \code{annotation}:the annotation of the region, it can be promoter,downstream,gene,UTR... 
#'              or the reason of liftOver FAILED.} }
#'
#'     
#' @return annotation data
#'
#' @seealso \code{\link{change_genomic_coordinate}}
#' 
"mouse_region"


#' @title a example output of \code{create_region_object }
#'
#' @description \code{human_obj <- create_region_object(human_met, human_region)}
#'
#' @format  A \code{list} object containing the two elements: \itemize{ \item{
#'   \code{met}: A list containing methylation region data, where each entry in
#'   the list is an \eqn{L_{i} X D} dimensional matrix, where \eqn{L_{i}}
#'   denotes the number of CpGs found in region \code{i}. The columns contain
#'   the following information: \enumerate{ \item{ 1st column: Contains the
#'   locations of CpGs relative to centre. Note that the actual locations are
#'   scaled to the [-1,1] region. } \item{ 2nd column: Contains the total 
#'   number of reads at each CpG location.} \item{ 3rd column: Contains the
#'   methylated reads at each CpG location.} }.
#'   Rownames of each matrix contain the actual CpG genomic
#'   coordinates as <chr>:<location>. } \item{ \code{anno}: The annotation
#'   object.} } Note: The lengths of \code{met} and \code{anno} should match.
#'
#'     
#' @return region_object
#'
#' @seealso \code{\link{create_region_object}}
#' 
"human_obj"


#' @title a example output of \code{create_region_object }
#'
#' @description \code{mouse_obj <- create_region_object(mouse_met, mouse_region)}
#'
#' @format  A \code{list} object containing the two elements: \itemize{ \item{
#'   \code{met}: A list containing methylation region data, where each entry in
#'   the list is an \eqn{L_{i} X D} dimensional matrix, where \eqn{L_{i}}
#'   denotes the number of CpGs found in region \code{i}. The columns contain
#'   the following information: \enumerate{ \item{ 1st column: Contains the
#'   locations of CpGs relative to centre. Note that the actual locations are
#'   scaled to the [-1,1] region. } \item{ 2nd column: Contains the total 
#'   number of reads at each CpG location.} \item{ 3rd column: Contains the
#'   methylated reads at each CpG location.} }.
#'   Rownames of each matrix contain the actual CpG genomic
#'   coordinates as <chr>:<location>. } \item{ \code{anno}: The annotation
#'   object.} } Note: The lengths of \code{met} and \code{anno} should match.
#'
#'     
#' @return region_object
#'
#' @seealso \code{\link{create_region_object}}
#' 
"mouse_obj"

#' @title a example out of \code{infer_profiles_vb }
#'
#' @description \code{human_basis_profile <- create_rbf_object(M = 8)
#' human_fit_profiles <- infer_profiles_vb(X = human_obj$met, model = "binomial",
#'                         basis = human_basis_profile, is_parallel = TRUE, vb_max_iter = 100)}
#'
#' @format An object of class \code{infer_profiles_vb_}"obs_model" with the
#'   following elements: \itemize{ \item{ \code{W}: An Nx(M+1) matrix with the
#'   optimized parameter values. Each row of the matrix corresponds to each
#'   element of the list X; if X is a matrix, then N = 1. The columns are of the
#'   same length as the parameter vector w (i.e. number of basis functions). }
#'   \item{ \code{W_Sigma}: A list with covariance matrices for each element row
#'   in W.} \item{ \code{basis}: The basis object. } \item{\code{nll_feat}: NLL
#'   fit feature.} \item{\code{rmse_feat}: RMSE fit feature.}
#'   \item{\code{coverage_feat}: CpG coverage feature.} \item{\code{lb_feat}:
#'   Lower Bound feature.}}
#'
#'     
#' @return infer_profiles_vb
#'
#' @seealso \code{\link{infer_profiles_vb}}
#' 
"human_fit_profiles"

#' @title a example out of \code{infer_profiles_vb }
#'
#' @description \code{human_basis_mean <- create_rbf_object(M = 0)
#' human_obj <- create_region_object(human_met, human_region)
#' human_fit_mean <- infer_profiles_vb(X = human_obj$met, model = "binomial",
#'                                     basis = human_basis_mean, is_parallel = TRUE, vb_max_iter = 100)}
#'
#' @format An object of class \code{infer_profiles_vb_}"obs_model" with the
#'   following elements: \itemize{ \item{ \code{W}: An Nx(M+1) matrix with the
#'   optimized parameter values. Each row of the matrix corresponds to each
#'   element of the list X; if X is a matrix, then N = 1. The columns are of the
#'   same length as the parameter vector w (i.e. number of basis functions). }
#'   \item{ \code{W_Sigma}: A list with covariance matrices for each element row
#'   in W.} \item{ \code{basis}: The basis object. } \item{\code{nll_feat}: NLL
#'   fit feature.} \item{\code{rmse_feat}: RMSE fit feature.}
#'   \item{\code{coverage_feat}: CpG coverage feature.} \item{\code{lb_feat}:
#'   Lower Bound feature.}}
#'
#'     
#' @return infer_profiles_vb
#'
#' @seealso \code{\link{infer_profiles_vb}}
#' 
#' 
"human_fit_mean"


#' @title a example out of \code{infer_profiles_vb }
#'
#' @description \code{mouse_basis_profile <- create_rbf_object(M = 8)
#' mouse_obj <- create_region_object(mouse_met, mouse_region)
#' mouse_fit_profiles <- infer_profiles_vb(X = mouse_obj$met, model = "binomial",
#'                                         basis = mouse_basis_profile, is_parallel = TRUE, vb_max_iter = 100)}
#'
#' @format An object of class \code{infer_profiles_vb_}"obs_model" with the
#'   following elements: \itemize{ \item{ \code{W}: An Nx(M+1) matrix with the
#'   optimized parameter values. Each row of the matrix corresponds to each
#'   element of the list X; if X is a matrix, then N = 1. The columns are of the
#'   same length as the parameter vector w (i.e. number of basis functions). }
#'   \item{ \code{W_Sigma}: A list with covariance matrices for each element row
#'   in W.} \item{ \code{basis}: The basis object. } \item{\code{nll_feat}: NLL
#'   fit feature.} \item{\code{rmse_feat}: RMSE fit feature.}
#'   \item{\code{coverage_feat}: CpG coverage feature.} \item{\code{lb_feat}:
#'   Lower Bound feature.}}
#'
#'     
#' @return infer_profiles_vb
#'
#' @seealso \code{\link{infer_profiles_vb}}
#' 
"mouse_fit_profiles"


#' @title a example out of \code{infer_profiles_vb }
#'
#' @description \code{mouse_basis_mean <- create_rbf_object(M = 0)
#' mouse_obj <- create_region_object(mouse_met, mouse_region)
#' mouse_fit_mean <- infer_profiles_vb(X = mouse_obj$met, model = "binomial",
#'                                     basis = mouse_basis_mean, is_parallel = TRUE, vb_max_iter = 100)}
#'
#' @format An object of class \code{infer_profiles_vb_}"obs_model" with the
#'   following elements: \itemize{ \item{ \code{W}: An Nx(M+1) matrix with the
#'   optimized parameter values. Each row of the matrix corresponds to each
#'   element of the list X; if X is a matrix, then N = 1. The columns are of the
#'   same length as the parameter vector w (i.e. number of basis functions). }
#'   \item{ \code{W_Sigma}: A list with covariance matrices for each element row
#'   in W.} \item{ \code{basis}: The basis object. } \item{\code{nll_feat}: NLL
#'   fit feature.} \item{\code{rmse_feat}: RMSE fit feature.}
#'   \item{\code{coverage_feat}: CpG coverage feature.} \item{\code{lb_feat}:
#'   Lower Bound feature.}}
#'
#'     
#' @return infer_profiles_vb
#'
#' @seealso \code{\link{infer_profiles_vb}}
#' 
"mouse_fit_mean"

#' @title a example output of \code{adjusted_cosine_similarity }
#'
#' @description \code{similarity <- adjusted_cosine_similarity(queryProfiles=human_fit_profiles,
#'                                                            subjectProfiles=mouse_fit_profiles)}
#'
#' @format A vector of Normalized adjusted cosine similarity [0-1], 
#'     
#' @return A vector of numeric
#'  Note: the more close to 1 indicates the profiles of region more similar,
#'  the more close to 0 indicates the profiles more different.
#'
#' @seealso \code{\link{adjusted_cosine_similarity}}
"similarity"