#' Plot estimated track, behavioural states and observations on a map. 
#' 
#' Takes a fitted \code{fitSSM} object and plots the observed (data) and estimated 
#' locations on a map. For the behavioural models (DCRWS, hDCRWS), the estimated
#' locations are coloured according to the posterior mean behavioural state estimates.
#' 
#' @param fit an output object from \code{fitSSM}
#' @param onemap  If FALSE (default) then each estimated track is plotted on a separate 
#' map, if TRUE then tracks are combined on a single map.
#' @return Observed locations are plotted as '+' symbols and estimated locations as filled
#' circles. Individual track id's (for DCRW and DCRWS models) are displayed at the top of 
#' each plot, but only when \code{onemap = FALSE}. The model specified in \code{fitSMM} is
#' also displayed at the top. Takes advantage of \code{ggplot2} plotting functions.
#' 
#' Currently, results from the hierarchical models (hDCRW, hDCRWS) can only be plotted on 
#' a combined map.
#' @examples
#' \dontrun{
#' data(ellie)
#' fit <- fitSSM(ellie, model = "DCRWS", tstep = 1, adapt = 2000, samples = 1000, 
#'               thin = 2, span = 0.1)
#' map_ssm(fit, onemap = TRUE)
#' 
#' fit.h <- fitSSM(ellie, model = "hDCRWS", tstep = 1, adapt = 2000, samples = 1000, 
#'                 thin = 2, span = 0.1)
#' map_ssm(fit.h)
#' }
#' @importFrom ggplot2 ggplot aes ggtitle geom_point scale_color_gradient2 xlab ylab aes_string
#' @importFrom ggplot2 fortify
#' @importFrom ggplot2 geom_polygon
#' @importFrom ggplot2 coord_cartesian
#' @export 

map_ssm <- function(fit, onemap = TRUE) 
{
  if(!is.null(fit$model)) map_ssm.h(fit, onemap)
  else {
    map_ssm.d(fit, onemap)
  }
}  

map_ssm.d <- function(fit, onemap) {
  countriesHigh <- NULL
  data(countriesHigh, package = "rworldxtra", envir = environment())
  wm <- fortify(countriesHigh)
  
  plt <- function(m) {
    xl <- extendrange(m$data$lon, f = 0.2)
    yl <- extendrange(m$data$lat, f = 0.2)
    
    p <- ggplot() + geom_polygon(data = wm, aes_string(x = "long", y = "lat", group = "group"), 
                                 fill = grey(0.3)) + 
      coord_cartesian(xlim = xl, ylim = yl) + xlab("Longitude") + ylab("Latitude")
    if(!onemap) {
      p <- p + ggtitle(paste(unique(as.character(m$summary$id)), "; ", m$model, sep = ""))
    }
    else if(onemap) {
      p <- p + ggtitle(m$model)
    }
    
    if(m$model == "DCRWS") {
      p <- p + geom_point(data = m$data, aes_string(x = "lon", y = "lat", group = NULL), 
                          colour = "goldenrod", pch = "+", size = 4) +
        geom_point(data = m$summary, aes_string(x = "lon", y = "lat", group = NULL, colour = "b"), size = 1.25) +
        scale_color_gradient2(midpoint = 1.5, low = "blue", mid = "white", high = "red")
    }
    
    else {
      p <- p + geom_point(data = m$data, aes_string(x = "lon", y = "lat", group = NULL), 
                          colour = "goldenrod", pch = "+", size = 4) +
        geom_point(data = m$summary, aes_string(x = "lon", y = "lat", group = NULL), colour = 'dodgerblue', 
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

map_ssm.h <- function(fit, onemap) {
  countriesHigh <- NULL
  data(countriesHigh, package = "rworldxtra", envir = environment())
  wm <- fortify(countriesHigh)
  
  xl <- extendrange(fit$data$lon, f = 0.2)
  yl <- extendrange(fit$data$lat, f = 0.2)
  
  p <- ggplot() + geom_polygon(data = wm, aes_string(x = "long", y = "lat", group = "group"), 
                               fill = grey(0.3)) + 
    coord_cartesian(xlim = xl, ylim = yl) + xlab("Longitude") + ylab("Latitude")
  if(!onemap) {
    p <- p + ggtitle(paste(unique(as.character(fit$summary$id)), "; ", fit$model, sep = ""))
  }
  else if(onemap) {
    p <- p + ggtitle(fit$model)
  }
 
  if(fit$model == "hDCRWS") {
    p <- p + geom_point(data = fit$data, aes_string(x = "lon", y = "lat", group = NULL), 
                        colour = "goldenrod", pch = "+", size = 4) +
      geom_point(data = fit$summary, aes_string(x = "lon", y = "lat", group = NULL, colour = "b"), size = 1.25) +
      scale_color_gradient2(midpoint = 1.5, low = "blue", mid = "white", high = "red")
  }
  
  else {
    p <- p + geom_point(data = fit$data, aes_string(x = "lon", y = "lat", group = NULL), 
                        colour = "goldenrod", pch = "+", size = 4) +
      geom_point(data = fit$summary, aes_string(x = "lon", y = "lat", group = NULL), colour = 'dodgerblue', 
                 size = 1.25)
  }
  p
}