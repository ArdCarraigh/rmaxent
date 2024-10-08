#' Calculate Multivariate Environmental Similarity
#'
#' Calculate Multivariate Environmental Similarity and most dissimilar/similar 
#' variables with respect to a reference dataset, for a set of environmental 
#' variables.
#'
#' @param x a `Raster*`, `list`, `matrix`, or `data.frame`
#'   where each layer/column/element represents focal values of an environmental
#'   variable.
#' @param ref a `list`, `matrix`, or `data.frame` where each
#'   column/element represents reference values for an environmental variable
#'   (corresponding to those given in `x`).
#' @param full (logical) should similarity values be returned for all variables?
#'   If `FALSE` (the default), then only the minimum similarity scores
#'   across variables will be returned.
#' @return If `x` is a `Raster*` object, this function returns a list 
#'   containing:
#'   - `similarity`: a `RasterStack` giving the environmental similarities for
#'   each variable in `x` (only included when `full=TRUE`); 
#'   - `similarity_min`: a `Raster` layer giving the minimum similarity value 
#'   across all variables for each location (i.e. the MESS);
#'   - `mod`: a factor `Raster` layer indicating which variable was most 
#'   dissimilar to its reference range (i.e. the MoD map, Elith et al. 2010); 
#'   and
#'   - `mos`: a factor `Raster` layer indicating which variable was most 
#'   similar to its reference range.
#'   
#'   If `x` is a `list`, `matrix`, or `data.frame`, the function will return
#'   a list as above, but with `RasterStack` and `Raster` objects replaced by 
#'   matrix and vectors.
#' @details `similarity` uses the MESS algorithm described in Appendix S3
#'   of Elith et al. 2010.
#' @keywords maxent, mess, similarity, environment
#' @references 
#' * Elith, J., Kearney, M., and Phillips, S. (2010) [The art of modelling range-shifting species](https://doi.org/10.1111/j.2041-210X.2010.00036.x). _Methods in Ecology and Evolution_, 1: 330-342. doi:10.1111/j.2041-210X.2010.00036.x
#' @importFrom raster stack nlayers init as.data.frame raster
#' @importFrom methods is
#' @importFrom stats na.omit
#' @export
#' @examples
#' library(dismo)
#' library(raster)
#' ff <- list.files(system.file('ex', package='dismo'), '\\.grd$', 
#'                  full.names=TRUE )
#' predictors <- stack(grep('biome', ff, value=TRUE, invert=TRUE))
#' occ <- read.csv(system.file('ex/bradypus.csv', package='dismo'))[, -1]
#' ref <- extract(predictors, occ)
#' mess <- similarity(predictors, ref, full=TRUE)
#' 
#' \dontrun{
#' library(rasterVis)
#' library(RColorBrewer)
#' levelplot(mess$mod, col.regions=brewer.pal(8, 'Set1'))
#' levelplot(mess$mos, col.regions=brewer.pal(8, 'Set1'))
#' }
similarity <- function (x, ref, full = FALSE) 
{
  if (!methods::is(ref, "data.frame")) {
    ref <- as.data.frame(ref)
  }
  if (is(x, "Raster")) {
    r <- TRUE
    if (isTRUE(full)) {
      out <- raster::stack(replicate(raster::nlayers(x), 
                                     raster::init(x, function(x) NA)))
    }
    else {
      out <- raster::init(x, function(x) NA)
    }
  }
  else r <- FALSE
  ref <- stats::na.omit(ref)
  
  if (!methods::is(x, "data.frame")) {
    factor_bool <- is.factor(x)
    x <- as.data.frame(values(x))
    x[,factor_bool] <- as.factor(x[,factor_bool])
  }
  
  new_order <- sapply(colnames(ref), function(name){which(name == colnames(x))})
  x <- x[,new_order]
  fact <- t(matrix(sapply(ref,is.factor)))
  
  if(any(!fact)){
    ref_numerical <- as.data.frame(ref[,!fact])
    colnames(ref_numerical) <- colnames(ref)[!fact]
    x_numerical <- as.data.frame(x[,!fact])
    colnames(x_numerical) <- colnames(x)[!fact]
    
    if (is.null(dim(ref_numerical))) {
      rng <- as.data.frame(range(ref_numerical, na.rm = TRUE))
    }
    else {
      rng <- as.data.frame(apply(ref_numerical, 2, range, na.rm = TRUE))
    }
    pct_less <- mapply(function(x, ref) {
      findInterval(x, sort(ref))/length(ref)
    }, x_numerical, ref_numerical, SIMPLIFY = FALSE)
    sim_numerical <- mapply(function(f, rng, p) {
      ifelse(f == 0, (p - rng[1])/diff(rng) * 100, 
             ifelse(f > 0 & f <= 0.5, f * 200, 
                    ifelse(f > 0.5 & f < 1, (1 - f) * 200, (rng[2] - p)/diff(rng) * 100)))
    }, pct_less, rng, x_numerical)
  }
  
  if(any(fact)){
    ref_categorical <- as.data.frame(as.character(ref[,fact]))
    colnames(ref_categorical) <- colnames(ref)[fact]
    x_categorical <- as.data.frame(as.character(x[,fact]))
    colnames(x_categorical) <- colnames(x)[fact]
    
    sim_categorical <- x_categorical
    for(i in 1:ncol(x_categorical)){
      temp_ref <- ref_categorical[,i]
      temp_x <- x_categorical[,i]
      temp_table <- table(temp_ref)
      sim_categorical[,i] <- sapply(temp_x, function(val){
        ifelse(is.na(val), NA, temp_table[val])})
      sim_categorical[,i] <- sim_categorical[,i]/length(temp_ref)
    }
    sim_categorical <- as.matrix(sim_categorical)
    sim_categorical[!is.na(sim_categorical) & sim_categorical == 0] <- -Inf
  }
  
  if(all(!fact)) sim <- sim_numerical
  else if(all(fact)) sim <- sim_categorical
  else sim <- cbind(sim_numerical, sim_categorical)
  
  min_sim <- if (is.matrix(sim)) 
    apply(sim, 1, min)
  else (min(sim))
  mins <- apply(sim, 1, which.min)
  most_dissimilar_vec <- unlist(ifelse(lengths(mins) == 0, NA, mins))
  maxs <- apply(sim, 1, which.max)
  most_similar_vec <- unlist(ifelse(lengths(maxs) == 0, NA, maxs))
  if (isTRUE(r)) {
    most_dissimilar <- raster::raster(out)
    most_dissimilar[] <- most_dissimilar_vec
    most_dissimilar <- as.factor(most_dissimilar)
    levels(most_dissimilar)[[1]] <- data.frame(ID = seq_len(ncol(sim)), 
                                               var = colnames(sim))
    most_similar <- raster::raster(out)
    most_similar[] <- most_similar_vec
    most_similar <- as.factor(most_similar)
    levels(most_similar)[[1]] <- data.frame(ID = seq_len(ncol(sim)), 
                                            var = colnames(sim))
    out_min <- raster::raster(out)
    out_min[] <- min_sim
    if (isTRUE(full)) {
      out[] <- sim
      list(similarity = out, similarity_min = out_min, 
           mod = most_dissimilar, mos = most_similar)
    }
    else list(similarity_min = out_min, mod = most_dissimilar, 
              mos = most_similar)
  }
  else {
    if (isTRUE(full)) {
      list(similarity = sim, similarity_min = min_sim, 
           mod = most_dissimilar_vec, mos = most_similar_vec)
    }
    else list(similarity_min = min_sim, mod = most_dissimilar_vec, 
              mos = most_similar_vec)
  }
}
