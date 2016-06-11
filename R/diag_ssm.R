#' Plot standard McMC convergence diagnostics to help determine lack of model convergence.
#' 
#' Takes a fitted \code{fit_ssm} object and uses standard McMC convergence diagnostic plots to
#' aid assessment of lack of convergence.
#' 
#' @param fit an output object from \code{fit_ssm}
#' @return Uses plotting functions from Martyn Plummer's \code{coda} package to help
#' diagnose lack of convergence for the core model parameters. The traceplot shows the time 
#' series for both McMC chains; the density plot shows the density estimate for each parameter;
#' the autocorrelation plots show the within-chain sample autocorrelation for each parameter;
#' the G-B-R shrink factor plot shows the evolution of Gelman and Rubin's shrink factor for 
#' increasing number of iterations. See the \code{coda} package for further details.
#' @references Brooks SP, Gelman A (1998) General methods for monitoring convergence of 
#' iterative simulations. Journal of Computational and Graphical Statistics 7:434-455
#' @examples
#' \dontrun{
#' data(ellie)
#' fit <- fit_ssm(ellie, model = "DCRWS", tstep = 1, adapt = 2000, samples = 1000, 
#'               thin = 2, span = 0.1)
#' diag_ssm(fit)
#' 
#' # increase burnin, posterior sample numbers, and thinning factor
#' fit2 <- fit_ssm(ellie, model = "DCRWS", tstep = 1, adapt = 5000, samples = 5000, 
#'               thin = 5, span = 0.1)
#' diag_ssm(fit2)
#' }             
#' @importFrom coda autocorr.plot traceplot gelman.plot as.mcmc.list nvar varnames<-
#' @export 
#' 
diag_ssm <- function(fit)
{
  if(!is.null(fit$model)) diag_ssm.h(fit)
  else {
    diag_ssm.d(fit)
  }
}  

diag_ssm.d <- function(fit) {
  
  dostuff <- function(m){
    md <- m$model
    if(md ==  "DCRW") layout(matrix(1:20, 4, 5, byrow = TRUE), widths = c(3, 3, 3, 3, 3), heights = rep(1, 5))
    else {
      layout(matrix(1:25, 5, 5, byrow = TRUE), widths = c(3, 3, 3, 3, 3), heights = rep(1, 5))
    }
    
    par(mar = c(1, 1, 1, 2), oma = c(3, 4, 2, 0))
    
    #gamma; single (if DCRW) or _1 and _2 (if DCRWS)
    foo <- as.mcmc.list(m$mcmc$gamma)
    varnames(foo) <- NULL
    if(nvar(foo) ==  2) {
      tmp <- list(as.mcmc.list(foo[, 1]), as.mcmc.list(foo[, 2]))
    }
    else {
      tmp <- list(foo)
    }
    t <- ifelse(md ==  "DCRWS", 2, 1)
    sapply(1:t, function(i) {
      #sample trace
      traceplot(tmp[[i]], col = c('blue', 'red'), lwd = 0.5, lty = c(1, 1), mgp = c(1, 0.5, 0), las = 1, 
                tcl = -0.2, cex = 0.7, smooth = TRUE, main = "")
      if(t > 1) mtext(bquote(gamma[.(i)]), 2, 3)
      else { mtext(expression(gamma), 2, 3) }
      if(i ==  1) mtext("Trace", 3, 0, cex = 0.7)	

      #posterior density
      foo1 = density(tmp[[i]][[1]], bw = 1.06 * min(sd(tmp[[i]][[1]]), IQR(tmp[[i]][[1]]) / 1.34) * 
                       length(tmp[[i]][[1]])^-0.2)
      foo1$x <- ifelse(foo1$x > 1, 1, foo1$x)
      foo1$x <- ifelse(foo1$x < 0, 0, foo1$x)
      foo2 = density(tmp[[i]][[2]], bw = 1.06 * min(sd(tmp[[i]][[2]]), IQR(tmp[[i]][[2]]) / 1.34) * 
                       length(tmp[[i]][[2]])^-0.2)
      foo2$x <- ifelse(foo2$x > 1, 1, foo2$x)
      foo2$x <- ifelse(foo2$x < 0, 0, foo2$x)
      plot(foo1, col = 'blue', lwd = 0.5, mgp = c(1, 0.5, 0), las = 1, tcl = -0.2, 
           cex = 0.7, ylim = c(-0.03 * max(foo1$y, foo2$y), max(foo1$y, foo2$y)), 
           xlim = c(min(foo1$x, foo2$x), max(foo1$x, foo2$x)), main = "")
      if(i ==  1) mtext("Density", 3, 0, cex = 0.7)	
      lines(foo2, col = 'red', lwd = 0.5)
      rug(tmp[[i]][[1]], lwd = 0.2, ticksize = 0.05, col = 'blue')
      rug(tmp[[i]][[2]], lwd = 0.2, ticksize = 0.05, col = 'red')	
      
      #sample autocorrelation
      autocorr.plot(tmp[[i]], auto.layout = FALSE, ask = FALSE, las = 1, tcl = -0.2, 
                    col = 'steelblue3', mgp = c(1, 0.2, 0))
      if(i ==  1) mtext("ACF: chain 1", 3, -1, outer = TRUE, adj = 0.5, cex = 0.7)
      if(i ==  1) mtext("ACF: chain 2", 3, -1, outer = TRUE, adj = 0.7, cex = 0.7)			
      if(i ==  1) mtext(paste(md, ": ", unique(m$summary$id), sep = ""), 3, 0.5, outer = TRUE)
      
      #G-B-R shrink factor (r-hat) plot
      gelman.plot(tmp[[i]], auto.layout = FALSE, ask = FALSE, xlab = "", ylab = "", 
                  las = 1, tcl = -0.2, mgp = c(1, 0.2, 0))
      if(i ==  1) mtext("G-R-B shrink factor", 3, -1, outer = TRUE, adj = 0.94, cex = 0.7)			
    })   
    
    #Sigma
    tmp = as.mcmc.list(m$mcmc$Sigma)
    varnames(tmp) <- NULL

    #longitude sigma 
    #sample trace
    traceplot(tmp[, 1], col = c('blue', 'red'), lwd = 0.5, lty = c(1, 1), mgp = c(1, 0.5, 0), las = 1, 
              tcl = -0.2, main = "", cex = 0.7, smooth = TRUE)
    mtext(expression(Sigma[lon]), 2, 3)	
    
    #posterior density
    foo1 = density(tmp[[1]][, 1], bw = 1.06 * min(sd(tmp[[1]][, 1]), IQR(tmp[[1]][, 1]) / 1.34)
                   * length(tmp[[1]][, 1])^-0.2)
    foo2 = density(tmp[[2]][, 1], bw = 1.06 * min(sd(tmp[[2]][, 1]), IQR(tmp[[2]][, 1]) / 1.34)
                   * length(tmp[[2]][, 1])^-0.2)			
    plot(foo1, col = 'blue', lwd = 0.5, mgp = c(1, 0.5, 0), las = 1, tcl = -0.2, main = "", 
         cex = 0.7, ylim = c(-0.03 * max(foo1$y, foo2$y), max(foo1$y, foo2$y)), 
         xlim = c(min(foo1$x, foo2$x), max(foo1$x, foo2$x)))
    lines(foo2, col = 'red', lwd = 0.5)
    rug(tmp[[1]][, 1], lwd = 0.2, ticksize = 0.05, col = 'blue')
    rug(tmp[[2]][, 1], lwd = 0.2, ticksize = 0.05, col = 'red')
    
    #sample autocorrelation
    autocorr.plot(tmp[, 1], auto.layout = FALSE, ask = FALSE, las = 1, tcl = -0.2, 
                  col = 'steelblue3', mgp = c(1, 0.2, 0))
    
    #G-B-R shrink factor (r-hat) plot
    gelman.plot(tmp[, 1], auto.layout = FALSE, ask = FALSE, xlab = "", ylab = "", 
                las = 1, tcl = -0.2, mgp = c(1, 0.2, 0))
    
    #latitude sigma
    #sample trace
    traceplot(tmp[, 4], col = c('blue', 'red'), lwd = 0.5, lty = c(1, 1), mgp = c(1, 0.5, 0), las = 1, 
              tcl = -0.2, main = "", cex = 0.7, smooth = TRUE)
    mtext(expression(Sigma[lat]), 2, 3)	
    
    #posterior density
    foo1 = density(tmp[[1]][, 4], bw = 1.06 * min(sd(tmp[[1]][, 4]), IQR(tmp[[1]][, 4]) / 1.34)
                   * length(tmp[[1]][, 1])^-0.2)
    foo2 = density(tmp[[2]][, 4], bw = 1.06 * min(sd(tmp[[2]][, 4]), IQR(tmp[[2]][, 4]) / 1.34)
                   * length(tmp[[2]][, 4])^-0.2)			
    plot(foo1, col = 'blue', lwd = 0.5, mgp = c(1, 0.5, 0), las = 1, tcl = -0.2, main = "", 
         cex = 0.7, ylim = c(-0.03 * max(foo1$y, foo2$y), max(foo1$y, foo2$y)), 
         xlim = c(min(foo1$x, foo2$x), max(foo1$x, foo2$x)))
    lines(foo2, col = 'red', lwd = 0.5)
    rug(tmp[[1]][, 4], lwd = 0.2, ticksize = 0.05, col = 'blue')
    rug(tmp[[2]][, 4], lwd = 0.2, ticksize = 0.05, col = 'red')	
    
    #sample autocorrelation
    autocorr.plot(tmp[, 4], auto.layout = FALSE, ask = FALSE, las = 1, tcl = -0.2, 
                  col = 'steelblue3', mgp = c(1, 0.2, 0))
    
    #G-B-R shrink factor (r-hat) plot
    gelman.plot(tmp[, 4], auto.layout = FALSE, ask = FALSE, xlab = "", ylab = "",  
                las = 1, tcl = -0.2, mgp = c(1, 0.2, 0))
    
    #psi
    tmp = as.mcmc.list(m$mcmc$psi)
    varnames(tmp) <- NULL
    
    #sample trace
    traceplot(tmp, col = c('blue', 'red'), lwd = 0.5, lty = c(1, 1), mgp = c(1, 0.5, 0), las = 1, 
              tcl = -0.2, main = "", cex = 0.7, smooth = TRUE)
    mtext(expression(psi), 2, 3)
    
    #posterior density
    foo1 = density(tmp[[1]], bw = 1.06 * min(sd(tmp[[1]][, 1]), IQR(tmp[[1]]) / 1.34) * 
                     length(tmp[[1]])^-0.2)
    foo2 = density(tmp[[2]], bw = 1.06 * min(sd(tmp[[2]][, 1]), IQR(tmp[[2]]) / 1.34) * 
                     length(tmp[[2]])^-0.2)			
    plot(foo1, col = 'blue', lwd = 0.5, mgp = c(1, 0.5, 0), las = 1, tcl = -0.2, main = "", 
         cex = 0.7, ylim = c(-0.03 * max(foo1$y, foo2$y), max(foo1$y, foo2$y)), 
         xlim = c(min(foo1$x, foo2$x), max(foo1$x, foo2$x)))
    lines(foo2, col = 'red', lwd = 0.5)	
    rug(tmp[[1]], lwd = 0.2, ticksize = 0.05, col = 'blue')
    rug(tmp[[2]], lwd = 0.2, ticksize = 0.05, col = 'red')
    
    #sample autocorrelation
    autocorr.plot(tmp, auto.layout = FALSE, ask = FALSE, las = 1, tcl = -0.2, 
                  col = 'steelblue3', mgp = c(1, 0.2, 0))
    
    #G-B-R shrink factor (r-hat) plot
    gelman.plot(tmp, auto.layout = FALSE, ask = FALSE, xlab = "", ylab = "", 
                las = 1, tcl = -0.2, mgp = c(1, 0.2, 0))
  }
  
  lapply(fit, dostuff)	
  invisible()
}

diag_ssm.h <- function(fit) {
  md = fit$model
  if(md ==  "hDCRW") layout(matrix(1:15, 3, 5, byrow = TRUE), widths = c(3, 3, 3, 3, 3), heights = rep(1, 4))
  else {
    layout(matrix(1:20, 4, 5, byrow = TRUE), widths = c(3, 3, 3, 3, 3), heights = rep(1, 5))
  }
  par(mar = c(1, 1, 1, 2), oma = c(3, 4, 2, 0))
  
  #gamma
  foo <- as.mcmc.list(fit$mcmc$gamma)
  varnames(foo) <- NULL
  if(nvar(foo) ==  2) {
    tmp <- list(as.mcmc.list(foo[, 1]), as.mcmc.list(foo[, 2]))
  }
  else {
    tmp <- list(foo)
  }
  t <- ifelse(md ==  "hDCRWS", 2, 1)
  sapply(1:t, function(i) {
    traceplot(tmp[[i]], col = c('blue', 'red'), lwd = 0.5, lty = c(1, 1), mgp = c(1, 0.5, 0), las = 1, 
              tcl = -0.2, cex = 0.7, smooth = TRUE, main = "")
    if(t > 1) mtext(bquote(gamma[.(i)]), 2, 3)
    else { mtext(expression(gamma), 2, 3) }
    if(i ==  1) mtext("Trace", 3, 0, cex = 0.7)	
    foo1 = density(tmp[[i]][[1]], bw = 1.06 * min(sd(tmp[[i]][[1]]), IQR(tmp[[i]][[1]]) / 1.34) * 
                     length(tmp[[i]][[1]])^-0.2)
    foo1$x <- ifelse(foo1$x > 1, 1, foo1$x)
    foo1$x <- ifelse(foo1$x < 0, 0, foo1$x)
    foo2 = density(tmp[[i]][[2]], bw = 1.06 * min(sd(tmp[[i]][[2]]), IQR(tmp[[i]][[2]]) / 1.34) * 
                     length(tmp[[i]][[2]])^-0.2)
    foo2$x <- ifelse(foo2$x > 1, 1, foo2$x)
    foo2$x <- ifelse(foo2$x < 0, 0, foo2$x)
    plot(foo1, col = 'blue', lwd = 0.5, mgp = c(1, 0.5, 0), las = 1, tcl = -0.2, 
         cex = 0.7, ylim = c(-0.03 * max(foo1$y, foo2$y), max(foo1$y, foo2$y)), 
         xlim = c(min(foo1$x, foo2$x), max(foo1$x, foo2$x)), main = "")
    if(i ==  1) mtext("Density", 3, 0, cex = 0.7)	
    lines(foo2, col = 'red', lwd = 0.5)
    rug(tmp[[i]][[1]], lwd = 0.2, ticksize = 0.05, col = 'blue')
    rug(tmp[[i]][[2]], lwd = 0.2, ticksize = 0.05, col = 'red')		
    autocorr.plot(tmp[[i]], auto.layout = FALSE, ask = FALSE, las = 1, tcl = -0.2, 
                  col = 'steelblue3', mgp = c(1, 0.2, 0))
    if(i ==  1) mtext("ACF: chain 1", 3, -1, outer = TRUE, adj = 0.5, cex = 0.7)
    if(i ==  1) mtext("ACF: chain 2", 3, -1, outer = TRUE, adj = 0.7, cex = 0.7)			
    if(i ==  1) mtext(md, 3, 0.5, outer = TRUE)
    
    #G-B-R shrink factor (r-hat) plot
    gelman.plot(tmp[[i]], auto.layout = FALSE, ask = FALSE, xlab = "", ylab = "", 
                las = 1, tcl = -0.2, mgp = c(1, 0.2, 0))
    if(i ==  1) mtext("G-R-B shrink factor", 3, -1, outer = TRUE, adj = 0.94, cex = 0.7)			
  })   
  
  #Sigma
  tmp = as.mcmc.list(fit$mcmc$Sigma)
  varnames(tmp) <- NULL
  
  #longitude sigma
  #sample trace
  traceplot(tmp[, 1], col = c('blue', 'red'), lwd = 0.5, lty = c(1, 1), mgp = c(1, 0.5, 0), las = 1, 
            tcl = -0.2, main = "", cex = 0.7, smooth = TRUE)
  mtext(expression(Sigma[lon]), 2, 3)	
  
  #posterior density
  foo1 = density(tmp[[1]][, 1], bw = 1.06 * min(sd(tmp[[1]][, 1]), IQR(tmp[[1]][, 1]) / 1.34)
                 * length(tmp[[1]][, 1])^-0.2)
  foo1$x <- ifelse(foo1$x < 0, 0, foo1$x)
  foo2 = density(tmp[[2]][, 1], bw = 1.06 * min(sd(tmp[[2]][, 1]), IQR(tmp[[2]][, 1]) / 1.34)
                 * length(tmp[[2]][, 1])^-0.2)
  foo2$x <- ifelse(foo2$x < 0, 0, foo2$x)
  plot(foo1, col = 'blue', lwd = 0.5, mgp = c(1, 0.5, 0), las = 1, tcl = -0.2, main = "", 
       cex = 0.7, ylim = c(-0.03 * max(foo1$y, foo2$y), max(foo1$y, foo2$y)), 
       xlim = c(min(foo1$x, foo2$x), max(foo1$x, foo2$x)))
  lines(foo2, col = 'red', lwd = 0.5)
  rug(tmp[[1]][, 1], lwd = 0.2, ticksize = 0.05, col = 'blue')
  rug(tmp[[2]][, 1], lwd = 0.2, ticksize = 0.05, col = 'red')	
  
  #sample autocorrelation
  autocorr.plot(tmp[, 1], auto.layout = FALSE, ask = FALSE, las = 1, tcl = -0.2, 
                col = 'steelblue3', mgp = c(1, 0.2, 0))
  
  #G-B-R shrink factor (r-hat) plot
  gelman.plot(tmp[, 1], auto.layout = FALSE, ask = FALSE, xlab = "", ylab = "", 
              las = 1, tcl = -0.2, mgp = c(1, 0.2, 0))
  
  #latitude sigma
  #sample trace
  traceplot(tmp[, 4], col = c('blue', 'red'), lwd = 0.5, lty = c(1, 1), mgp = c(1, 0.5, 0), las = 1, 
            tcl = -0.2, main = "", cex = 0.7, smooth = TRUE)
  mtext(expression(Sigma[lat]), 2, 3)	
  
  #posterior density
  foo1 = density(tmp[[1]][, 4], bw = 1.06 * min(sd(tmp[[1]][, 4]), IQR(tmp[[1]][, 4]) / 1.34)
                 * length(tmp[[1]][, 1])^-0.2)
  foo1$x <- ifelse(foo1$x < 0, 0, foo1$x)
  foo2 = density(tmp[[2]][, 4], bw = 1.06 * min(sd(tmp[[2]][, 4]), IQR(tmp[[2]][, 4]) / 1.34)
                 * length(tmp[[2]][, 4])^-0.2)
  foo2$x <- ifelse(foo2$x < 0, 0, foo2$x)
  plot(foo1, col = 'blue', lwd = 0.5, mgp = c(1, 0.5, 0), las = 1, tcl = -0.2, main = "", 
       cex = 0.7, ylim = c(-0.03 * max(foo1$y, foo2$y), max(foo1$y, foo2$y)), 
       xlim = c(min(foo1$x, foo2$x), max(foo1$x, foo2$x)))
  lines(foo2, col = 'red', lwd = 0.5)
  rug(tmp[[1]][, 4], lwd = 0.2, ticksize = 0.05, col = 'blue')
  rug(tmp[[2]][, 4], lwd = 0.2, ticksize = 0.05, col = 'red')	
  
  #sample autocorrelation
  autocorr.plot(tmp[, 4], auto.layout = FALSE, ask = FALSE, las = 1, tcl = -0.2, 
                col = 'steelblue3', mgp = c(1, 0.2, 0))
  
  #G-B-R shrink factor (r-hat) plot
  gelman.plot(tmp[, 4], auto.layout = FALSE, ask = FALSE, xlab = "", ylab = "",  
              las = 1, tcl = -0.2, mgp = c(1, 0.2, 0))
  
  #psi
  #one psi for each individual dataset (N total datasets)
  N <- fit$N
  foo <- as.mcmc.list(fit$mcmc$psi)
  varnames(foo) <- NULL
  if(N > 1) {
    tmp <- lapply(1:N, function(i){
      as.mcmc.list(foo[, i])
    })
  }
  else {
    tmp <- list(foo)
  }
  ids <- as.character(unique(fit$summary$id))
  layout(matrix(1:(N * 5), N, 5, byrow = TRUE), widths = c(3, 3, 3, 3, 3), heights = rep(1, N + 1))
  sapply(1:N, function(i) {
    
    # sample trace
    traceplot(tmp[[i]], col = c('blue', 'red'), lwd = 0.5, lty = c(1, 1), mgp = c(1, 0.5, 0), las = 1, 
              tcl = -0.2, main = "", cex = 0.7, smooth = TRUE)
    mtext(bquote(psi[.(ids[i])]), 2, 3)
    if(i ==  1) mtext("Trace", 3, 0, cex = 0.7)
    
    #posterior density
    foo1 <- density(tmp[[i]][[1]], bw = 1.06 * min(sd(tmp[[i]][[1]]), IQR(tmp[[i]][[1]]) / 1.34) * 
                     length(tmp[[i]][[1]])^-0.2)
    foo1$x <- ifelse(foo1$x < 0, 0, foo1$x)
    foo2 <- density(tmp[[i]][[2]], bw = 1.06 * min(sd(tmp[[i]][[2]]), IQR(tmp[[i]][[2]]) / 1.34) * 
                     length(tmp[[i]][[2]])^-0.2)
    foo2$x <- ifelse(foo2$x < 0, 0, foo2$x)
    plot(foo1, col = 'blue', lwd = 0.5, mgp = c(1, 0.5, 0), las = 1, tcl = -0.2, main = "", 
         cex = 0.7, ylim = c(-0.03 * max(foo1$y, foo2$y), max(foo1$y, foo2$y)), 
         xlim = c(min(foo1$x, foo2$x), max(foo1$x, foo2$x)))
    lines(foo2, col = 'red', lwd = 0.5)	
    rug(tmp[[i]][[1]], lwd = 0.2, ticksize = 0.05, col = 'blue')
    rug(tmp[[i]][[2]], lwd = 0.2, ticksize = 0.05, col = 'red')	
    if(i ==  1) mtext("Density", 3, 0, cex = 0.7)
    
    #sample autocorrelation
    autocorr.plot(tmp[[i]], auto.layout = FALSE, ask = FALSE, las = 1, tcl = -0.2, 
                  col = 'steelblue3', mgp = c(1, 0.2, 0)) 
    if(i ==  1) mtext("ACF: chain 1", 3, -1, outer = TRUE, adj = 0.5, cex = 0.7)
    if(i ==  1) mtext("ACF: chain 2", 3, -1, outer = TRUE, adj = 0.7, cex = 0.7)	
    
    #G-B-R shrink factor (r-hat) plot
    gelman.plot(tmp[[i]], auto.layout = FALSE, ask = FALSE, xlab = "", ylab = "", 
                las = 1, tcl = -0.2, mgp = c(1, 0.2, 0))
    if(i ==  1) mtext("G-R-B shrink factor", 3, -1, outer = TRUE, adj = 0.94, cex = 0.7)
    if(i ==  1) mtext(md, 3, 0.5, outer = TRUE)
  })
  invisible()
}