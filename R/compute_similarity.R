#' @title compute the similarity of queryProfiles and subjectProfiles
#' @description compute the similarity of queryProfiles and subjectProfiles by adjusted cosine 
#'              distance and return a vector of numeric match the queryProfiles and subjectProfiles.
#' @param queryProfiles a result of \code{infer_profiles_vb}
#' @param subjectProfiles a result of \code{infer_profiles_vb}
#' @param is_parallel Logical, indicating if code should be run in parallel.
#' @param no_cores Number of cores given as input
#' 
#' @return a vector of Normalized adjusted cosine similarity [0-1], 
#'  Note: the more close to 1 indicates the profiles of region more similar,
#'  the more close to 0 indicates the profiles more different.
#' 
#' @author Hongen Kang  \email{geneprophet@163.com}
#' @export
#'
#' @examples
#' \dontrun{
#' similarity <- adjusted_cosine_similarity(queryProfiles=human_fit_profiles, 
#'                                          subjectProfiles=mouse_fit_profiles)
#'                                          
#' #check the number of non-NA
#' length(which(!is.na(similarity)))
#' 
#' #check the number of >0.9
#' length(which(similarity>0.9))
#' 
#' #check the number of <0.1
#' length(which(similarity<0.1))
#' }
#' 
adjusted_cosine_similarity <- function(queryProfiles,
                                       subjectProfiles,
                                       is_parallel = TRUE,
                                       no_cores = NULL){
  
  assertthat::assert_that(methods::is(queryProfiles,"infer_profiles"))
  assertthat::assert_that(methods::is(subjectProfiles,"infer_profiles"))
  if(NROW(queryProfiles$W)!=NROW(subjectProfiles$W)){
    stop("ERROR:the queryProfiles and subjectProfiles must have the same length!")
  }
  assertthat::assert_that(NROW(queryProfiles$W)==NROW(subjectProfiles$W))
  N = NROW(queryProfiles$W)
  if(is_parallel){
    n_cores <- .parallel_cores(is_parallel = is_parallel,no_cores = no_cores)
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    result <- foreach::`%dopar%`(obj = foreach::foreach(i=1:N,.packages = "BSDMR",.inorder = TRUE,
                                                        .combine = 'c',.multicombine = TRUE, .maxcombine = 100),
                                 ex = {out <- .compute_similarity(queryProfiles = queryProfiles,subjectProfiles = subjectProfiles,index = i)})
    # Stop parallel execution
    parallel::stopCluster(cl)
    doParallel::stopImplicitCluster()
  }else{
    result <- foreach::`%do%`(obj = foreach::foreach(i=1:N,.packages = "BSDMR",.inorder = TRUE,
                                                     .combine = 'c',.multicombine = TRUE, .maxcombine = 100),
                              ex = {out <- .compute_similarity(queryProfiles = queryProfiles,subjectProfiles = subjectProfiles,index = i)})
  }
  
 
  # 
  # Result <- matrix(data = result,ncol = 5,byrow = T)
  # final_result <- data.frame(similarity=Result[,1],queryProfilesIndex = Result[,2],subjectProfilesIndex=Result[,3],sim = Result[,4],corr = Result[,5])
  
  return(result)
  
}

