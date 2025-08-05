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
similarity <- function (x, ref, full = FALSE, filename='', ...) 
{
  if (!methods::is(ref, "data.frame")) {
    ref <- as.data.frame(ref)
  }
  ref <- stats::na.omit(ref)
  
  if (methods::is(x, "Raster")) {
    r <- TRUE
    out <- raster(x)
    filename <- trim(filename)
    nms <- names(x)
    x <- as.data.frame(raster::values(x))
    colnames(x) <- nms
  }
  else r <- FALSE
  
  if (!methods::is(x, "data.frame")) {
    x <- as.data.frame(x)
  }
  
  x <- x[colnames(x) %in% colnames(ref)]
  ref <- ref[colnames(ref) %in% colnames(x)]
  ref <- ref[colnames(x)]
  fact <- c(t(matrix(sapply(ref,is.factor))))
  nms <- colnames(x)
  nl <- ncol(x)
  
  ref_numerical <- as.data.frame(ref[!fact])
  rng <- as.data.frame(apply(ref_numerical, 2, range, na.rm = TRUE))
  ref_categorical <- as.data.frame(ref[fact])
  
  .mess <- function(x, ref_numerical, rng, ref_categorical, fact){
    if(any(!fact)){
      x_numerical <- as.data.frame(x[!fact])
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
      x_categorical <- as.data.frame(x[fact])
      sim_categorical <- x_categorical
      for(i in 1:ncol(x_categorical)){
        temp_ref <- ref_categorical[,i]
        temp_x <- x_categorical[,i]
        temp_table <- table(temp_ref)
        sim_categorical[,i] <- sapply(temp_x, function(val){
          ifelse(is.na(val), NA, ifelse(val %in% names(temp_table), temp_table[val], -Inf))})
        sim_categorical[,i] <- sim_categorical[,i]/length(temp_ref) * 100
      }
      sim_categorical <- as.matrix(sim_categorical)
    }
    
    if(all(!fact)) sim <- as.data.frame(sim_numerical)
    else if(all(fact)) sim <- as.data.frame(sim_categorical)
    else sim <- as.data.frame(cbind(sim_numerical, sim_categorical))
    sim <- sim[colnames(ref)]
    sim[!is.na(sim) & sim == -Inf]  <- min(0, min(sim[!is.na(sim) & sim != -Inf]) - 1)
    min_sim <- apply(sim, 1, min)
    mins <- apply(sim, 1, which.min)
    most_dissimilar_vec <- unlist(ifelse(lengths(mins) == 0, NA, mins))
    maxs <- apply(sim, 1, which.max)
    most_similar_vec <- unlist(ifelse(lengths(maxs) == 0, NA, maxs))
    list(similarity = sim, similarity_min = min_sim, 
         mod = most_dissimilar_vec, mos = most_similar_vec)
  }
  
  if (isTRUE(r)) {
    out <- brick(out, nl=nl+3)
    if(filename == ''){
      vv <- .mess(x, ref_numerical, rng, ref_categorical, fact)
      out <- setValues(out, as.matrix(cbind(vv$similarity, vv$similarity_min, unname(vv$mod), unname(vv$mos))))
    }
    else{
      tr <- blockSize(out)
      pb <- pbCreate(tr$n, ...)	
      out <- writeStart(out, filename, ...)
      for (i in 1:tr$n) {
        start <- ncol(out) * (tr$row[i]-1) + 1
        end <- start + (ncol(out) * tr$nrows[i]-1)
        vv <- x[start:end,]
        vv <- .mess(vv, ref_numerical, rng, ref_categorical, fact)
        out <- writeValues(out, as.matrix(cbind(vv$similarity, vv$similarity_min, unname(vv$mod), unname(vv$mos))), tr$row[i])
        pbStep(pb) 
      }
      out <- writeStop(out)
      pbClose(pb)
    }
    names(out) <- c(nms, "similarity_min", "mod", "mos")
    most_dissimilar <- as.factor(out[['mod']])
    levels(most_dissimilar)[[1]] <- data.frame(ID = seq_len(nl), var = nms)
    most_similar <- as.factor(out[['mos']])
    levels(most_similar)[[1]] <- data.frame(ID = seq_len(nl), var = nms)
    out_min <- out[['similarity_min']]
    out <- out[[nms]]
    
    if (isTRUE(full)) {
      out <- list(similarity = out, similarity_min = out_min, 
                  mod = most_dissimilar, mos = most_similar)
    }
    else out <- list(similarity_min = out_min, mod = most_dissimilar, 
                     mos = most_similar)
  }
  
  else {
    out <- .mess(x, ref_numerical, rng, ref_categorical, fact)
    out$mod <- as.factor(out$mod)
    levels(out$mod) <- nms
    out$mos <- as.factor(out$mos)
    levels(out$mos) <- nms
    if (!isTRUE(full)) {
      out <- out[c('similarity_min', 'mod', 'mos')]
    }
  }
  out
}
