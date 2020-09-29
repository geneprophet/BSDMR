#' @name create_region_object
#' 
#' @aliases create_region_obj
#'
#' @title Create genomic region data
#'
#' @description \code{create_region_object} creates genomic regions (e.g. forms
#'   methylation regions data) using as input methylation and annotation data
#'   with genomic regions of interest.
#'
#' @param met_dt A \code{GRanges} object with methylation data, whose format
#'   should be similar to \code{\link{read_methylation_report}}.
#' @param anno_dt A \code{GRanges} object with annotation data, whose format
#'   should be similar to  \code{\link{read_annotation}}.
#' @param cov Integer defining the minimum coverage of CpGs that each region
#'   must contain.
#' @param ignore_strand Logical, whether or not to ignore strand information.
#' @param is_parallel Logical, indicating if code should be run in parallel.
#' @param no_cores if you have specify the is_parallel is ture ,you can specify the number of parallel cores
#'
#' @return A \code{list} object containing the two elements: \itemize{ \item{
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
#' @author Hongen Kang  \email{geneprophet@163.com}
#'
#' @examples
#' \dontrun{
#' # Download the files and change the working directory to that location
#' human_met <- read_methylation_report("name_of_met_file")
#' human_anno <- read_annotation("name_of_anno_file")
#' human_obj <- create_region_object(human_met, human_anno)
#'}
#' @seealso \code{\link{read_methylation_report}}, \code{\link{read_annotation}}
#'
#' @importFrom methods is
#' @export
create_region_object <- function(met_dt, 
                                 anno_dt,
                                 cov = 5, 
                                 ignore_strand = TRUE, 
                                 is_parallel = TRUE,
                                 no_cores = NULL) {
  fmin = -1
  fmax = 1
  sd_thresh = 1e-1000
  message("Creating methylation regions ...")
  assertthat::assert_that(methods::is(met_dt, "GRanges"))
  assertthat::assert_that(methods::is(anno_dt, "GRanges"))
  # Find overlaps between met and anno
  overlaps <- GenomicRanges::findOverlaps(query = anno_dt, subject = met_dt, ignore.strand = T)
  query_hits <- S4Vectors::queryHits(overlaps)
  subj_hits  <- S4Vectors::subjectHits(overlaps)
  gen_ind    <- unique(query_hits)
  centre     <- anno_dt$center     # Central locations
  id         <- anno_dt$id         # (Ensembl) IDs
  chrom      <- as.character(anno_dt@seqnames) # Chrom info
  strand     <- as.character(GenomicRanges::strand(anno_dt))
  cpg_loc    <- GenomicRanges::ranges(met_dt)@start  # CpG locations
  met <- met_dt$methylated_reads    # Methylated read
  
  # Number of columns for each matrix
  if (is.null(met_dt$total_reads)) {
    D <- 2
  } else {
    total <- met_dt$total_reads
    D <- 3
  }
  
  # Extract upstream and downstream lengths
  N <- NROW(anno_dt)
  up_anno <-
    vector(mode = "integer", N)    # Start location in chromosome
  down_anno <-
    vector(mode = "integer", N)  # End location in chromosome
  anno_start <- GenomicRanges::ranges(anno_dt)@start # Start info
  anno_end <-
    anno_start + GenomicRanges::ranges(anno_dt)@width - 1 # End info
  for (i in 1:N) {
    # Depending on the strand we change regions up or downstream of centre
    if (identical(strand[i], "-")) {
      up_anno[i] <- centre[i] - anno_end[i]
      down_anno[i] <- abs(anno_start[i] - centre[i])
    } else {
      up_anno[i] <- anno_start[i] - centre[i]
      down_anno[i] <- anno_end[i] - centre[i]
    }
  }
  rm(anno_start, anno_end)
  
  if(is_parallel){
    n_cores <- .parallel_cores(is_parallel = is_parallel,no_cores = no_cores)
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    result <- foreach::`%dopar%`(obj = foreach::foreach(i=1:N,.inorder = TRUE,.packages = "BSDMR",
                                                        .multicombine = TRUE, .maxcombine = 1000),
                                 ex = {out <- .create_met_region(met_dt=met_dt,gen_ind=gen_ind,query_hits = query_hits,subj_hits = subj_hits, cov = cov,D=D,met=met,total=total, sd_thresh=sd_thresh,cpg_loc=cpg_loc,centre=centre,chrom=chrom,strand=strand,up_anno=up_anno,down_anno=down_anno,fmin=fmin,fmax=fmax, index = i)})
    # Stop parallel execution
    parallel::stopCluster(cl)
    doParallel::stopImplicitCluster()
  } else {
    result <- foreach::`%do%`(obj = foreach::foreach(i=1:N,.packages = "BSDMR",.inorder = TRUE,
                                                     .multicombine = TRUE, .maxcombine = 1000),
                              ex = {out <- .create_met_region(met_dt=met_dt,gen_ind=gen_ind,query_hits = query_hits,subj_hits = subj_hits, cov = cov,D=D,met=met,total=total, sd_thresh=sd_thresh,cpg_loc=cpg_loc,centre=centre,chrom=chrom,strand=strand,up_anno=up_anno,down_anno=down_anno,fmin=fmin,fmax=fmax,index = i)})
  }
  
  met_region = result


  return(structure(list(met = met_region, anno = anno_dt),
                   class = "region_object"))
}


#' @name create_basis
#' 
#' @aliases basis create_basis_function
#'
#' @title Create basis objects
#'
#' @description These functions create different basis objects, which can be
#'   used as input to complex functions in order to perform computations
#'   depending on the class of the basis function.
#'
#' @param M The number of the basis functions. In case of Fourier basis, this
#'   number should be even, since we need to have pairs of sines and cosines and
#'   the constant term is added by default.
#' @param gamma Inverse width of radial basis function.
#' @param mus Optional centers of the RBF.
#' @param eq_spaced_mus Logical, if TRUE, equally spaced centers are created,
#'   otherwise centers are created using \code{\link[stats]{kmeans}} algorithm.
#' @param whole_region Logical, indicating if the centers will be evaluated
#'   equally spaced on the whole region, or between the min and max of the
#'   observation values.
#'
#' @return A basis object of class 'rbf'.
#'
#' @author Hongen Kang  \email{geneprophet@163.com}
#'
#' @seealso  \code{\link{design_matrix}}
#' 
#'
#' @examples
#' human_basis_profile <- create_rbf_object(M = 8)
#' human_basis_mean <- create_rbf_object(M = 0)
#' #---------------------------------
#'
#' @export
create_rbf_object <- function(M = 8, 
                              gamma = NULL,
                              mus = NULL,
                              eq_spaced_mus = TRUE, 
                              whole_region = TRUE){
  # Check that M is numberic and integer
  assertthat::assert_that(is.numeric(M))
  assertthat::assert_that(is.logical(eq_spaced_mus))
  assertthat::assert_that(is.logical(whole_region))
  assertthat::assert_that(M %% 1 == 0)
  assertthat::assert_that(M > -1)
  if (!is.null(gamma)) {
    assertthat::assert_that(is.numeric(gamma))
    assertthat::assert_that(gamma > 0)
  }else {gamma <- M ^ 2 / (abs(1) + abs(-1)) ^ 2 }
  if (!is.null(mus)) {
    assertthat::assert_that(is.vector(mus))
    assertthat::assert_that(M == length(mus))
  }else{
    if (eq_spaced_mus) {
      mus <- vector(mode = "numeric", M)
      if (whole_region) {
        for (i in 1:M) { mus[i] <- i * ((1 - (-1)) / (M + 1) ) + (-1) }
      }
    }
  }
  obj <- structure(list(M = M, mus = mus, gamma = gamma,
                        eq_spaced_mus = eq_spaced_mus,
                        whole_region = whole_region),
                   class = "rbf")
  return(obj)
}


