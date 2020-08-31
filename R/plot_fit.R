#' Plot the 1-D time-series of estimated location and behavioural states
#' 
#' Takes a fitted \code{fit_ssm} object and plots the observed (data), estimated 
#' location and behavioural states (posterior means) as 1-D time-series. Each
#' individual dataset is plotted separately.
#' 
#' @param fit an output object from \code{fit_ssm}
#' @return Observed locations are plotted as filled circles and estimated locations as blue
#' lines with the 95\% credible interval as a ribbon. Uses \code{ggplot2} plotting functions.
#' 
#' @examples
#' \dontrun{
#' data(ellie)
#' fit.s <- fit_ssm(ellie, model = "DCRWS", tstep = 1, adapt = 100, samples = 100, 
#'               thin = 1, span = 0.1)
#' plot_fit(fit.s)
#' 
#' hfit.s <- fit_ssm(ellie, model = "hDCRWS", tstep = 1, adapt = 100, samples = 100, 
#'                 thin = 1, span = 0.1)
#' plot_fit(hfit.s)
#' }
#' @importFrom ggplot2 ggplot aes ggtitle geom_point scale_color_gradient2 xlab ylab aes_string 
#' @importFrom ggplot2 ylim geom_line geom_ribbon
#' @importFrom gridExtra grid.arrange
#' @export 

plot_fit <- function(fit) 
{

  plt <- function(d) {
    #longitude
    yl <- range(c(d$data$lon, d$summary$lon.025, d$summary$lon.975), na.rm = TRUE)
    p1 <- ggplot() + geom_point(data = d$data, aes_string(x = "date", y = "lon", group = NULL), 
                                 colour = "firebrick", size = 0.75) + 
      ylab("Longitude") + xlab("") + ylim(yl[1], yl[2]) + 
      geom_line(data = d$summary, aes_string(x = "date", y = "lon", group = NULL), 
                          colour = "dodgerblue") + 
      geom_ribbon(data = d$summary, aes_string(x = "date", ymin = "lon.025", ymax = "lon.975"), 
                  fill = "dodgerblue", alpha = 0.5) + 
      ggtitle(paste(unique(as.character(d$summary$id)), "; ", d$model, sep = ""))
    
    #latitude
    yl <- range(c(d$data$lat, d$summary$lat.025, d$summary$lat.975), na.rm = TRUE)
    p2 <- ggplot() + geom_point(data = d$data, aes_string(x = "date", y = "lat", group = NULL), 
                                colour = "firebrick", size = 0.75) + 
      ylab("Latitude") + xlab("") + ylim(yl[1], yl[2]) +
      geom_line(data = d$summary, aes_string(x = "date", y = "lat", group = NULL), 
                            colour = "dodgerblue") + 
      geom_ribbon(data = d$summary, aes_string(x = "date", ymin = "lat.025", ymax = "lat.975"), 
                  fill = "dodgerblue", alpha = 0.5)
    
    #behaviour
    if(d$model == "DCRWS" || d$model == "hDCRWS") {
      p3 <- ggplot(data = d$summary) + geom_line(aes_string(x = "date", y = "b", group = NULL), 
                                                 colour = "dodgerblue", size = 1) +
        ylab("Behavioural state") + xlab("")
    }
    
  if(d$model == "DCRW" || d$model == "hDCRW") {  
    grid.arrange(p1, p2, heights = c(2, 2))
  }
    else {
      grid.arrange(p1, p2, p3, heights = c(2, 2, 2))
    }
  }  
  
  if(class(fit) == "ssm") lapply(fit, plt)
  
  else if(class(fit) == "hssm") {
    N <- fit$N
    s <- with(fit, split(summary, summary$id))
    d <- with(fit, split(data, data$id))
    fit <- lapply(1:N, function(i) {
      list(summary = s[[i]], data = d[[i]], model = fit$model)
    })
    lapply(fit, plt)
  }
  
  invisible()  
}