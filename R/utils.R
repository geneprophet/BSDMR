# (Internal) Apply the radial basis function
#
# Apply the RBF function to the input X.
#
# @param X Input data.
# @param mus Centers from where we should compute the distance of the data X.
# @param gamma Inverse width of radial basis function.
#
# @return Input X, after being transformed from the RBF.
#
# @author Hongen Kang  \email{geneprophet@163.com}
#
#
.rbf_basis <- function(X,
                       mus, 
                       gamma = 1){
  return(exp( (-1) * gamma * sum( (X - mus) ^ 2) ))
}


#' @name design_matrix
#' @rdname design_matrix
#' @aliases designmatrix des_matrix des_mat
#'
#' @title Generic function for creating design matrices
#'
#' @description These functions call the appropriate methods depending on the
#'   class of the object \code{obj} to create RBF design matrices.
#'
#' @param obj A basis function object.
#' @param obs A vector of observations.
#' @param ... Additional parameters.
#'
#' @return A design matrix object
#'
#' @author Hongen Kang  \email{geneprophet@163.com}
#'
#' @export
design_matrix <- function(obj, 
                          ...){
  UseMethod("design_matrix")
}


#' @rdname design_matrix
#'
design_matrix.default <- function(obj, 
                                  ...){
  stop(paste0("Object type '", class(obj), "' is not implemented."))
}

#' @rdname design_matrix
#'
#' @export
design_matrix.rbf <- function(obj,
                              obs,
                              ...){
  assertthat::assert_that(methods::is(obj, "rbf"))
  
  if (is.na(obs[1])) {
    return(list(H = NA, basis = obj))
  } else {
    if(methods::is(obs,"matrix")){
      obs = obs[,1]
    }
    assertthat::assert_that(is.vector(obs))
    N   <- length(obs)  # Length of the dataset
    if (obj$M == 0) { H <- matrix(1, nrow = N, ncol = 1); obj$mus <- 0;
    }else{
      if (is.null(obj$mus)) {
        if (obj$eq_spaced_mus) {
          obj$mus <- vector(mode = "numeric", obj$M)
          if (!obj$whole_region) {
            # TODO: Keep this as functionality?
            for (i in 1:obj$M) {
              obj$mus[i] <- i*(max(obs) - min(obs))/(obj$M + 1) +
                min(obs)
            }
          }
        }else {
          repeat {
            # TODO: Keep this as functionality?
            km <- stats::kmeans(obs, obj$M, iter.max = 30, nstart = 10)
            if (min(km$size) > 0) { break } # Accept non-empty clusters
          }
          obj$mus <- km$centers  # RBF centers
        }
      }
      # Convert the 'obs' vector to an N x 1 dimensional matrix
      obs <- as.matrix(obs)
      H <- matrix(1, nrow = N, ncol = obj$M + 1)
      for (j in 1:obj$M) {
        H[, j + 1] <- apply(obs,1,.rbf_basis,mus = obj$mus[j],
                            gamma = obj$gamma)
      }
    }
    return(list(H = H, basis = obj))
  }
}



# Center CpG locations relative to TSS
#
# \code{center_loc} centera CpG locations relative to TSS
#
# @param region CpG locations
# @param tss TSS location
# @param strand_direction Strand direction
#
# @return Centered location data relative to TSS
#
.do_centre_loc <- function(region,
                           centre, 
                           strand_direction){
  assertthat::assert_that(is.character(strand_direction))
  center <- region - centre
  if (identical(strand_direction, "-")){
    center  <- (-center)  # If '-' strand, swap CpG locations
  }
  return(center)
}
# scaling_location
# Compute the min-max scaling
#
# \code{.minmax_scaling} normalizes a given vector using the the min-max
# scaling method. More formally:
# \deqn{scaled = \frac{data -x_{min}}{x_{max} - x_{min}} \times (f_{max} -
#  f_{min}) + f_{min}}
#
# @param data Vector with numeric data to be scaled.
# @param xmin Optional minimum value, otherwise \code{min(data)} will be used.
# @param xmax Optional maximum value, otherwise \code{max(data)} will be used.
# @param fmin Optional minimum range value, default is -1.
# @param fmax Optional maximum range value, default is 1.
#
# @return The scaled data in the given range, default is between (-1, 1). If
#  xmin = xmax the input vector \code{data} is returned.
#
.minmax_scaling <- function(data,
                            xmin = NULL, 
                            xmax = NULL,
                            fmin = -1, 
                            fmax = 1){
  if (is.null(xmin)) { xmin <- min(data) }
  if (is.null(xmax)) { xmax <- max(data) }
  if ( (xmin - xmax) == 0) { return(data) }
  minmax <- (data - xmin) / (xmax - xmin)
  minmax_scaled <- minmax * (fmax - fmin) + fmin
  return(minmax_scaled)
}


# @title Number of parallel cores
#
# @description Function for creating the number of parallel cores that will be
#   used during EM.
# @param no_cores Number of cores given as input
# @param is_parallel Logical, did we require parallel computations
# @param M Total number of sources
#
.parallel_cores <- function(no_cores=NULL,
                            is_parallel=FALSE, 
                            max_cores = NULL){
  if (is_parallel) { # If parallel mode is ON
    # If number of cores is not given
    if (is.null(no_cores)) { no_cores <- parallel::detectCores() - 1
    } else{if (no_cores >= parallel::detectCores()) {
      no_cores <- parallel::detectCores() - 1 }
    }
    if (is.na(no_cores)) { no_cores <- 2 }
    if (!is.null(max_cores)) {
      if (no_cores > max_cores) { no_cores <- max_cores }
    }
  }
  return(no_cores)
}
# Compute predictive distribution of inferred profiles using VB
.predictive_infer_profile <- function(model,
                                      x,
                                      region = 1){
  # Predictive mean
  H <- design_matrix(model$basis, x)$H
  if(methods::is(model$W[region, ], "numeric")){
    W_pred <- c(H %*% model$W[region, ])
  }else{
    W_pred <- c(H %*% as.vector(model$W[region, ]$w1))
  }
  if (methods::is(model, "infer_profiles_vb_binomial") ||
      methods::is(model, "infer_profiles_vb_bernoulli")) {
    # Predictive variance
    W_sd_pred <- sqrt(1 + diag(H %*% model$W_Sigma[[region]] %*% t(H)))
    W_pred <- pnorm(W_pred / W_sd_pred)
  }else {stop("Not predictive distribution for this object") }
  return(list(W_pred = W_pred, W_sd_pred = W_sd_pred))
}
# @title findOverlap
# 根据hits聚点成region,间隔小于200bp,连续10个
# parallel function 
.parallel_find_region <- function(methylation,
                                  annotation,
                                  index,
                                  min_sites_number,
                                  max_distance_between_sites){
  anno <- annotation[index]
  hits <- GenomicRanges::findOverlaps(methylation,anno,ignore.strand=TRUE)
  begin <- vector(mode = "double", length = 0)
  end <- vector(mode = "double", length = 0)
  anno_index <- vector(mode = "double", length = 0)
  if(length(hits)>=min_sites_number){
    count = 1
    b = GenomicRanges::start(methylation[S4Vectors::queryHits(hits[1])])
    e = NULL
    for (i in 2:length(hits)) {
      #current <- GenomicRanges::start(methylation[S4Vectors::queryHits(hits[i])])
      #before <- GenomicRanges::start(methylation[S4Vectors::queryHits(hits[i-1])])
      if((GenomicRanges::start(methylation[S4Vectors::queryHits(hits[i])])-GenomicRanges::start(methylation[S4Vectors::queryHits(hits[i-1])]))<max_distance_between_sites){
        count = count + 1
        print(count)
        #finsh the iteration when reach the end of CoG_gr
        if(count==length(hits)){
          e = GenomicRanges::start(methylation[S4Vectors::queryHits(hits[i])])
          begin <- c(begin,b)
          end <- c(end,e)
          #anno_index <- c(anno_index,index)
        }
      }else{
        if(count >= min_sites_number){
          e = GenomicRanges::start(methylation[S4Vectors::queryHits(hits[i-1])])
          begin <- c(begin,b)
          end <- c(end,e)
          #anno_index <- c(anno_index,index)
          b = GenomicRanges::start(methylation[S4Vectors::queryHits(hits[i])])
          count = 1
          e = NULL
        }else{
          count = 1
          b = GenomicRanges::start(methylation[S4Vectors::queryHits(hits[i])])
          e = NULL
        }
      }
    }
    
  }

  #if no region meet the suitation then return NULL
  if(length(end)==0){
    return(NULL)
  }else{
    anno_index <- rep(index,length(end))
    region <- data.frame(start=begin,end=end,anno_index=anno_index)
    return(region)
  }
 
}

.constructGRanges <- function(inputRegion,
                              out,
                              index){
  #判断Deleted
  if(length(out[[index]])==0){
    result <- GenomicRanges::GRanges(
      seqnames = "FAILED",
      ranges = IRanges::IRanges(start = 0,end = 0),
      strand = "-",
      id =  "liftOver FAILED",
      center = NA,
      annotation = "Deleted in new: Sequence intersects no chains"
    )
    return(result)
  }else{
    #判断split
    if(length(BiocGenerics::unique(BiocGenerics::as.vector(GenomicRanges::seqnames(out[[index]]))))>1 || length(BiocGenerics::unique(BiocGenerics::as.vector(BiocGenerics::strand(out[[index]]))))>1){
      result <- GenomicRanges::GRanges(
        seqnames = "FAILED",
        ranges = IRanges::IRanges(start = 0,end = 0),
        strand = "-",
        id =  "liftOver FAILED",
        center = NA,
        annotation = "Split in new: Sequence insufficiently intersects multiple chains"
      )
      return(result)
    }else{
      if(tail(BiocGenerics::as.vector(GenomicRanges::end(out[[index]])),1) > head(BiocGenerics::as.vector(GenomicRanges::start(out[[index]])),1)){
        result <- GenomicRanges::GRanges(
          seqnames = BiocGenerics::unique(BiocGenerics::as.vector(GenomicRanges::seqnames(out[[index]]))),
          ranges = IRanges::IRanges(start = head(BiocGenerics::as.vector(GenomicRanges::start(out[[index]])),1), end = tail(BiocGenerics::as.vector(GenomicRanges::end(out[[index]])),1)),
          strand = BiocGenerics::unique(BiocGenerics::as.vector(BiocGenerics::strand(out[[index]]))),
          id = "liftOver SUCCESS",
          center = round((head(BiocGenerics::as.vector(GenomicRanges::start(out[[index]])),1) + tail(BiocGenerics::as.vector(GenomicRanges::end(out[[index]])),1))/2),
          annotation = "Successed"
        )
        #判断Partially deleted
        if((GenomicRanges::width(inputRegion[index])-GenomicRanges::width(result))>30){
          result <- GenomicRanges::GRanges(
                  seqnames = "FAILED",
                  ranges = IRanges::IRanges(start = 0,end = 0),
                  strand = "-",
                  id = "liftOver FAILED",
                  center = NA,
                  annotation = "Partially deleted in new: Sequence insufficiently intersects one chain"
                )
        }
        return(result)
      }else{
        result <- GenomicRanges::GRanges(
          seqnames = BiocGenerics::unique(BiocGenerics::as.vector(GenomicRanges::seqnames(out[[index]]))),
          ranges = IRanges::IRanges(start =  tail(BiocGenerics::as.vector(GenomicRanges::end(out[[index]])),1), end = head(BiocGenerics::as.vector(GenomicRanges::start(out[[index]])),1)),
          strand = BiocGenerics::unique(BiocGenerics::as.vector(BiocGenerics::strand(out[[index]]))),
          id = "liftOver SUCCESS",
          center = round((tail(BiocGenerics::as.vector(GenomicRanges::end(out[[index]])),1) + head(BiocGenerics::as.vector(GenomicRanges::start(out[[index]])),1))/2),
          annotation = "Successed"
        )
        if((GenomicRanges::width(inputRegion[index])-GenomicRanges::width(result))>30){
          result <- GenomicRanges::GRanges(
            seqnames = "FAILED",
            ranges = IRanges::IRanges(start = 0,end = 0),
            strand = "-",
            id = "liftOver FAILED",
            center = NA,
            annotation = "Partially deleted in new: Sequence insufficiently intersects one chain"
          )
        }
        return(result)
      }
    }
  }
}


.create_met_region <- function(met_dt,
                               gen_ind,
                               query_hits,
                               subj_hits,
                               cov,
                               D,
                               met,
                               total,
                               sd_thresh,
                               cpg_loc,
                               centre,
                               chrom,
                               strand,
                               up_anno,
                               down_anno,
                               fmin,
                               fmax,
                               index){
  
  cpg_sites = subj_hits[which(query_hits==index)]
  
  if(length(cpg_sites) > cov){
    # If sd of the methylation level is above threshold
    if (D == 2) {
      obs_var <- stats::sd(met[cpg_sites])
    } else {
      obs_var <- stats::sd(met[cpg_sites] / total[cpg_sites])
    }
    if (obs_var > sd_thresh) {
      # Locations of CpGs in the genome
      region <- cpg_loc[cpg_sites]
      # Middle location for region[index]
      middle <- centre[gen_ind[which(gen_ind==index)]]
      # Extract chromosome information
      chrom_info <- chrom[gen_ind[which(gen_ind==index)]]
      # Extract strand information, i.e. direction
      strand_direction <- strand[gen_ind[which(gen_ind==index)]]
      # Extract upstream information
      upstream <- up_anno[gen_ind[which(gen_ind==index)]]
      # Extract downstream information
      downstream <- down_anno[gen_ind[which(gen_ind==index)]]
      # Shift CpG locations relative to TSS
      center_data  <- .do_centre_loc(
        region = region,
        centre = middle,
        strand_direction = strand_direction
      )
      # In "-" strand the order of the locations should change
      Order <- base::order(center_data)
      res <- matrix(data = 0, nrow = length(cpg_sites), ncol = D)
      # Store actual genomic coordinates of CpGs
      rownames(res) <- paste(chrom_info,
                                           region[Order], sep = ":")
      # Store normalized locations of CpGs
      res[, 1] <- round(
        .minmax_scaling(
          data = center_data[Order],
          xmin = upstream,
          xmax = downstream,
          fmin = fmin,
          fmax = fmax
        ),
        10
      )
      # Store reads in the corresponding locations
      if (D == 2) {
        res[, 2] <- met[cpg_sites][Order]
      } else{
        # Store total reads in the corresponding locations
        res[, 2] <- total[cpg_sites][Order]
        # Store methylated reads in the corresponding locations
        res[, 3] <- met[cpg_sites][Order]
      }
    }
    return(res)
  }else{
    res = NA
    return(res)
  }
}


.compute_similarity <- function(queryProfiles,
                                subjectProfiles,
                                index){
  if( is.na(queryProfiles$W[index,]) || is.na(subjectProfiles$W[index,])){
    return(NA)
  }else {
    a <- queryProfiles$W[index,]
    b <- subjectProfiles$W[index,]
    if(methods::is(a, "list")){
      a=as.vector(a$w1)
    }
    if(methods::is(b, "list")){
      b=as.vector(b$w1)
    }
    
    similarity <- crossprod(a-mean(a),b-mean(b))/sqrt(crossprod(a-mean(a))*crossprod(b-mean(b)))
    #sim <- crossprod(a,b)/sqrt(crossprod(a)*crossprod(b))
    #c <- cor(a,b)
    #normalize to [0-1]
    s <- 0.5 + 0.5*similarity
    #res <- c(s, index, which(subjectObj$anno$correspondingIndex == index),sim,c)
    return(s)
  }
}
