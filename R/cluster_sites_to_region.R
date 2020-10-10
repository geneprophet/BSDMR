#' @title cluster sites to region
#' 
#' @description cluster sites which meet the required suitations to region according to the annotation file
#' 
#' @param methylation \code{GRanges} object created by \code{read_methylation_report}.
#' @param annotation \code{GRanges} object created by \code{read_annotation}.
#' @param min_sites_number the minmal number of C sites in a region. the default is 10
#' @param ignore_strand Logical, whether or not to ignore strand information.
#' @param max_distance_between_sites the maxmal distance between C sites within a region.
#'                                   the unit of distance is unit and the default is 200
#' @param is_parallel Logical, indicating if code should be run in parallel.
#' @param no_cores if you have specify the is_parallel is ture ,you can specify the number of parallel cores
#' 
#' @return A \code{GRanges} object.
#' 
#' @author Hongen Kang  \email{geneprophet@163.com}
#' @export
#' @seealso \code{\link{read_methylation_report}},  \code{\link{read_annotation}}
#' @examples
#' \dontrun{
#' human_region <- cluster_sites_to_region(methylation = human_met,
#'                                         annotation = human_anno,is_parallel = TRUE)
#' }
#' 
#' 
cluster_sites_to_region <- function(methylation,
                                    annotation,
                                    min_sites_number = 10,
                                    max_distance_between_sites = 200,
                                    ignore_strand = TRUE,
                                    is_parallel = TRUE,
                                    no_cores = NULL){
  n_cores <- .parallel_cores(is_parallel = is_parallel,no_cores = no_cores)
  if(is_parallel){
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    # Parallel cluster sites to regions for each annotation region i.
    res <- foreach::`%dopar%`(obj = foreach::foreach(i=1:length(annotation),.packages = "BSDMR",.inorder = TRUE,
                                                     .combine = 'rbind',.multicombine = TRUE, .maxcombine = 100),
                              ex = {out <- .parallel_find_region(methylation = methylation,annotation = annotation,index = i,min_sites_number = min_sites_number,max_distance_between_sites = max_distance_between_sites,ignore_strand = ignore_strand)})
    # Stop parallel execution
    parallel::stopCluster(cl)
    doParallel::stopImplicitCluster()
  }else{
    res <- foreach::`%do%`(obj = foreach::foreach(i=1:length(annotation), .packages = "BSDMR",.inorder = TRUE,
                                                  .combine = 'rbind',.multicombine = TRUE, .maxcombine = 100),
                           ex = {out <- .parallel_find_region(methylation = methylation,annotation = annotation,index = i,min_sites_number = min_sites_number,max_distance_between_sites = max_distance_between_sites,ignore_strand = ignore_strand)})
  }
  
  result <- GenomicRanges::GRanges(
    seqnames = BiocGenerics::as.vector(GenomicRanges::seqnames(annotation[res$anno_index])),
    ranges = IRanges::IRanges(start = res$start, end = res$end),
    strand = BiocGenerics::as.vector(GenomicRanges::strand(annotation[res$anno_index])),
    id =  BiocGenerics::as.vector(annotation[res$anno_index]$id),
    center = round((res$start+res$end)/2),
    annotation = BiocGenerics::as.vector(annotation[res$anno_index]$annotation)
  )
  
  return(result)
}
