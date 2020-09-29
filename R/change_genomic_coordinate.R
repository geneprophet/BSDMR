#' @title change the genomic coordinate system of inputRegion by liftOver with chainFile
#'
#' @param inputRegion the region you want to change the genomic coordinate system,a \code{GRanges} object.
#'                    Generally, the inputRegin is a result of \code{cluster_sites_to_region}.
#' @param chainFile a chainFile of liftOver. for example, if you want mapping the the inputRegion from hg38 to mm10,
#'                 then the chainFile must be "hg38ToMm10.over.chain".
#'                  you can download from the source
#'                   <http://hgdownload.soe.ucsc.edu/downloads.html#liftover>
#' @param targetAnnotation the \code{GRanges} object to annotate the result of liftOver, a \code{GRanges} object.
#'                         Generally, the targetAnnotation is a result of \code{read_annotation}.
#' @param is_parallel Logical, indicating if code should be run in parallel.
#' @param no_cores if you have specify the is_parallel is ture ,you can specify the number of parallel cores
#' 
#' @return A \code{GRanges} object corrsponding to the inputRegion in another genomic coordinate system
#' 
#' @author Hongen Kang  \email{geneprophet@163.com}
#' @export
#'
#' @examples  
#' \dontrun{
#' filePath <- system.file("extdata", "hg38ToMm10.over.chain", package = "BSDMR")
#' mouse_region <- change_genomic_coordinate(human_region,filePath,mouse_anno)
#' }
#' 
change_genomic_coordinate <- function(inputRegion,
                                      chainFile,
                                      targetAnnotation,
                                      is_parallel = TRUE,
                                      no_cores = NULL){
  if(!methods::is(inputRegion,"GRanges")){
    stop(paste0(inputRegion,"ERROR:the paramater of inputRegion is not GRanges"))
  }
  
  chain <- rtracklayer::import.chain(chainFile)
  gr <- inputRegion
  over <- rtracklayer::liftOver(gr,chain)
  n_cores <- .parallel_cores(is_parallel = is_parallel,no_cores = no_cores)
  if(length(inputRegion)!=length(over)){
    stop("ERROR:the length of liftOver result list is not equal to the length of inputRegion")
  }
  
  suppressWarnings({
    if(is_parallel){
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)
      result <- foreach::`%dopar%`(obj = foreach::foreach(i=1:length(inputRegion),.packages = "BSDMR",.inorder = TRUE,
                                                          .combine = 'c',.multicombine = TRUE, .maxcombine = 100),
                                   ex = {out <- .constructGRanges(inputRegion = inputRegion,out = over,index = i)})
      # Stop parallel execution
      parallel::stopCluster(cl)
      doParallel::stopImplicitCluster()
    }else {
      result <- foreach::`%do%`(obj = foreach::foreach(i=1:length(inputRegion),  .packages = "BSDMR",.inorder = TRUE,
                                                       .combine = "c",.multicombine = TRUE, .maxcombine = 100),
                                ex = {out <- .constructGRanges(inputRegion = inputRegion,out = over,index = i)})
    }
    
    #把转换后的region注释到对应物种的基因组上
    hits <- GenomicRanges::findOverlaps(result,targetAnnotation,ignore.strand=T)
    
    #把转换成功的region注释到targetAnnotation
    for (j in 1:length(hits)) {
      result[S4Vectors::queryHits(hits[j])]$id = targetAnnotation[S4Vectors::subjectHits(hits[j])]$id
      result[S4Vectors::queryHits(hits[j])]$annotation = targetAnnotation[S4Vectors::subjectHits(hits[j])]$annotation
    }
    
  })
  
  return(result)
}

