#' Plot estimated track, behavioural states and observations on a map. 
#' 
#' Estimated location states are shaded by estimated behavioural state 
#' (for DCRWS, hDCRWS output). Observations are added as '+' symbols. 
#' For multi-individual output objects there is an option to include all
#' on a single map or to plot each individually.
#' 
#' @param fit an output object from \code{fitSSM}
#' @param onemap  logical (default is FALSE) indicating if individual tracks
#' are combined on a single map (TRUE) or on separate maps (FALSE)
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 fortify
#' @importFrom ggplot2 geom_polygon
#' @importFrom ggplot2 coord_cartesian
#' @export 

tplot <- function(fit, onemap = TRUE) 
{
  data(countriesHigh, package="rworldxtra")
  wm <- fortify(countriesHigh)

  plt <- function(m) {
    xl <- extendrange(m$data$lon, f = 0.2)
    yl <- extendrange(m$data$lat, f = 0.2)
    
    p <- ggplot() + geom_polygon(data = wm, aes(x = long, y = lat, group = group), 
                                 fill = grey(0.3)) + 
      coord_cartesian(xlim = xl, ylim = yl) 
    if(!onemap) {
      p <- p + ggtitle(paste(unique(as.character(m$summary$id)), "; ", m$model, sep = ""))
    }
    else if(onemap) {
      p <- p + ggtitle(m$model)
    }
    
    if(m$model == "DCRWS") {
      p <- p + geom_point(data = m$data, aes(x = lon, y = lat, group = NULL), 
                          colour = "goldenrod", pch = "+", size = 4) +
        geom_point(data = m$summary, aes(x = lon, y = lat, group = NULL, colour = b), size = 1.25) +
        scale_color_gradient2(midpoint = 1.5, low = "blue", mid = "white", high = "red")
    }
    
    else {
      p <- p + geom_point(data = m$data, aes(x = lon, y = lat, group = NULL), 
                          colour = "goldenrod", pch = "+", size = 4) +
        geom_point(data = m$summary, aes(x = lon, y = lat, group = NULL), colour = 'dodgerblue', 
                   size = 1.25)
    }
    p
  }  
if(length(fit) == 1) plt(fit[[1]])
  
else if(onemap && length(fit) > 1) {
  s <- do.call(rbind, lapply(fit, function(x) x$summary))
  d <- do.call(rbind, lapply(fit, function(x) x$data))
  m <- sapply(fit, function(x) x$model)[1]
  fit <- list(summary = s, data = d, model = m)
  plt(fit)
}  
  
else if(!onemap) {
  lapply(fit, plt)
}
  
}