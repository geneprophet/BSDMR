# Define ggplot2 theme
.gg_theme <- function(){
  p <- theme(
    plot.title = element_text(size = 18,face = 'bold',
                              margin = margin(0,0,3,0), hjust = 0.5),
    axis.text = element_text(size = rel(1.25), color = 'black'),
    axis.title = element_text(size = rel(1.45), color = 'black'),
    axis.title.y = element_text(margin = margin(0,10,0,0)),
    axis.title.x = element_text(margin = margin(10,0,0,0)),
    axis.ticks.x = element_line(colour = "black", size = rel(0.8)),
    axis.ticks.y = element_blank(),
    legend.position = "right",
    legend.key.size = unit(1.2, 'lines'),
    legend.title = element_text(size = 12, face = 'bold'),
    legend.text = element_text(size = 10),
    panel.border = element_blank(),
    panel.grid.major = element_line(colour = "gainsboro"),
    panel.background = element_blank()
  )
  return(p)
}

#' @title Plot inferred methylation profiles across a region
#'
#' @description Function for plotting the inferred methylation profiles across a
#'   given region, and optionally the mean methylation rate together with the
#'   observed methylation data, using \code{\link{ggplot2}}.
#'
#' @param region Genomic region number
#' @param obj_prof Inferred profile, i.e. output from
#'   \code{\link{infer_profiles_vb}}
#' @param obj_mean Inferred mean function, i.e. output from
#'   \code{\link{infer_profiles_vb}}
#' @param obs a list result of \code{\link{create_region_object}}
#' @param title Plot title
#' @param ... Additional parameters
#'
#' @return A ggplot2 object.
#'
#' @author Hongen Kang  \email{geneprophet@163.com}
#'
#' @examples
#' \dontrun{
#' # Fit methylation profiles using 8 RBFs and 0 RBF
#' human_basis_profile <- create_rbf_object(M = 8)
#' human_basis_mean <- create_rbf_object(M = 0)
#' human_obj <- create_region_object(human_met, human_region)
#' human_fit_profiles <- infer_profiles_vb(X = human_obj$met, model = "binomial",
#'    basis = human_basis_profile, is_parallel = TRUE, vb_max_iter = 100)
#' human_fit_mean <- infer_profiles_vb(X = human_obj$met, model = "binomial",
#'    basis = human_basis_mean, is_parallel = TRUE, vb_max_iter = 100)
#' p <- plot_infer_profiles(region = 44, obj_prof = human_fit_profiles,obj_mean = human_fit_mean, 
#'    obs = human_obj, title = paste0("Gene ID ",human_obj$anno$id[44]))
#' }
#' 
#' @seealso \code{\link{infer_profiles_vb}}, \code{\link{create_rbf_object}}, \code{\link{create_region_object}}
#' 
#' @export
#' 
plot_infer_profiles <- function(region = 1, obj_prof, obj_mean = NULL,
                                obs = NULL, title = "Inferred profiles",
                                ...) {
  if (is.na(obj_prof$W[region])){
    stop("Error:No methylation profile for this region!")
  }
  aes_xs <- seq(from = -1, to = 1, by = 0.001)
  x_axis = as.vector(GenomicRanges::seqnames(obs$anno[region]))
  y_axis = "methylation level"
  #@param x_labels x axis ticks labels
  #x_labels = c("Upstream", "", "Centre", "", "Downstream")
  x_labels = c(GenomicRanges::start(obs$anno[region]),"",obs$anno[region]$center,"",GenomicRanges::end(obs$anno[region]))
  obs = obs$met
  # For RMD CHECK to pass without NOTEs
  aes_ys = x = y = Model = ys_prof = ys_mean <- NULL
  ys_low = ys_high <- 0
  ys_low_prof = ys_high_prof = ys_low_mean = ys_high_mean <- 0
  if (methods::is(obj_prof, "infer_profiles_vb")) {
    tmp <- .predictive_infer_profile(obj_prof, aes_xs, region)
    ys_prof <- tmp$W_pred
    if (methods::is(obj_prof, "infer_profiles_vb_binomial") ||
        methods::is(obj_prof, "infer_profiles_vb_bernoulli")) {
      ys_low_prof <- ys_prof - ys_prof*(1 - ys_prof);
      ys_high_prof <- ys_prof + ys_prof*(1 - ys_prof)
    }
    
    if (!is.null(obj_mean)) {
      tmp <- .predictive_infer_profile(obj_mean, aes_xs, region)
      ys_mean <- tmp$W_pred
      if (methods::is(obj_mean, "infer_profiles_vb_binomial") ||
          methods::is(obj_mean, "infer_profiles_vb_bernoulli")) {
        ys_low_mean <- ys_mean - ys_mean*(1 - ys_mean);
        ys_high_mean <- ys_mean + ys_mean*(1 - ys_mean)
      }
    }
  }else{
    stop("No plotting function for this model!")
  }
  
  dt <- data.table::data.table(aes_xs = aes_xs, aes_ys = ys_prof,
                               ys_low = ys_low_prof, ys_high = ys_high_prof)
  if (is.null(obj_mean)) {
    p <- ggplot(dt, aes(x = aes_xs, y = aes_ys)) +
      geom_line(aes(x = aes_xs, y = aes_ys), size = 1.5, col = "darkblue")
    if (methods::is(obj_prof, "infer_profiles_vb") ||
        methods::is(obj_prof, "infer_profiles_gibbs") ) {
      p <- p + geom_ribbon(dt, mapping = aes(ymin = ys_low,
                                             ymax = ys_high), alpha = 0.2, size = 0.1, fill = "cornflowerblue")
    }
  }else{
    dt <- dt %>% .[, c("Model") := list("Profile")]
    dt_mean <- data.table::data.table(aes_xs = aes_xs, aes_ys = ys_mean,
                                      ys_low = ys_low_mean,
                                      ys_high = ys_high_mean,
                                      Model = "Mean")
    dt <- rbind(dt, dt_mean)
    if (methods::is(obj_prof, "infer_profiles_vb") ||
        methods::is(obj_prof, "infer_profiles_gibbs") ) {
      p <- ggplot(dt, aes(x = aes_xs, y = aes_ys, color = Model)) +
        geom_line(size = 1.2) +
        geom_ribbon(dt, mapping = aes(ymin = ys_low, ymax = ys_high,
                                      fill = Model),
                    alpha = 0.2, size = 0.1) +
        scale_color_brewer(palette = "Dark2") +
        scale_fill_brewer(palette = "Dark2")
    }else{
      p <- ggplot(dt, aes(x = aes_xs, y = aes_ys, color = Model)) +
        geom_line(size = 1.5) +
        scale_color_brewer(palette = "Dark2") +
        scale_fill_brewer(palette = "Dark2")
    }
  }
  if (!is.null(obs)) {
    # If we are given a list object
    if (is.list(obs)) { obs <- obs[[region]] }
    # If we have binomial observations
    if (NCOL(obs) == 3) {
      dt_obs <- data.table::data.table(x = obs[,1], y = obs[, 3]/obs[, 2])
    }else{
      dt_obs <- data.table::data.table(x = obs[,1], y = obs[, 2])
    }
    p <- p + geom_point(data = dt_obs, mapping = aes(x = x, y = y),
                        shape = 1, color = "blue", size = 2)
  }
  p <- p + scale_x_continuous(limits = c(-1, 1), labels = x_labels) +
    labs(title = title, x = x_axis, y = y_axis) +
    .gg_theme()
  return(p)
}

