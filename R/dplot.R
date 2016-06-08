#' dplot
#' 
#' @param fit fit
#'
#' @importFrom coda autocorr.plot traceplot as.mcmc.list
#' @export 
dplot = function(fit){

  dostuff = function(m){
    md = m$model
    
    layout(matrix(1:16,4,4,byrow=TRUE), widths=c(3,3,3,3), heights=rep(1,5))
    par(mar=c(1,1,1,2), oma=c(3,4,2,0))
    
    #gamma
    tmp = as.mcmc.list(m$mcmc$gamma)
    traceplot(tmp, col=c('blue','red'), lwd=0.5, lty=c(1,1), mgp=c(1,0.5,0), las=1,
              tcl=-0.2, cex=0.7, smooth=TRUE, main = "")
    mtext(expression(gamma), 2, 3)	
    mtext("Trace", 3, 0, cex=0.7)	
    foo1 = density(tmp[[1]], bw = 1.06 * min(sd(tmp[[1]][,1]), IQR(tmp[[1]])/1.34) * 
                     length(tmp[[1]])^-0.2)
    foo2 = density(tmp[[2]], bw = 1.06 * min(sd(tmp[[2]][,1]), IQR(tmp[[2]])/1.34) * 
                     length(tmp[[2]])^-0.2)			
    plot(foo1, col='blue', lwd=0.5, mgp=c(1,0.5,0), las=1, tcl=-0.2,
         cex=0.7, ylim=c(-0.03*max(foo1$y,foo2$y), max(foo1$y,foo2$y)), 
         xlim=c(min(foo1$x,foo2$x),max(foo1$x,foo2$x)), main="")
    mtext("Density", 3, 0, cex=0.7)	
    lines(foo2, col='red', lwd=0.5)
    rug(tmp[[1]], lwd=0.2, ticksize=0.05, col='blue')
    rug(tmp[[2]], lwd=0.2, ticksize=0.05, col='red')							
    autocorr.plot(tmp, auto.layout=FALSE, ask=FALSE, las=1, tcl=-0.2,
                  col='steelblue3', mgp=c(1,0.5,0), main = "")
    mtext("ACF: chain 1", 3, -1, outer=T, adj=0.635, cex=0.7)
    mtext("ACF: chain 2", 3, -1, outer=T, adj=0.91, cex=0.7)			
    mtext(paste(md, ": ", unique(m$summary$id), sep=""), 3, 0.5, outer=TRUE)
    
    #Sigma
    tmp = as.mcmc.list(m$mcmc$Sigma)
    traceplot(tmp[,1], col=c('blue','red'), lwd=0.5, lty=c(1,1), mgp=c(1,0.5,0), las=1,
              tcl=-0.2, main="", cex=0.7, smooth=TRUE)
    mtext(expression(Sigma[lon]), 2, 3)	
    foo1 = density(tmp[[1]][,1], bw = 1.06 * min(sd(tmp[[1]][,1]), IQR(tmp[[1]][,1])/1.34)
                   * length(tmp[[1]][,1])^-0.2)
    foo2 = density(tmp[[2]][,1], bw = 1.06 * min(sd(tmp[[2]][,1]), IQR(tmp[[2]][,1])/1.34)
                   * length(tmp[[2]][,1])^-0.2)			
    plot(foo1, col='blue', lwd=0.5, mgp=c(1,0.5,0), las=1, tcl=-0.2, main="",
         cex=0.7, ylim=c(-0.03*max(foo1$y,foo2$y), max(foo1$y,foo2$y)), 
         xlim=c(min(foo1$x,foo2$x),max(foo1$x,foo2$x)))
    lines(foo2, col='red', lwd=0.5)
    rug(tmp[[1]][,1], lwd=0.2, ticksize=0.05, col='blue')
    rug(tmp[[2]][,1], lwd=0.2, ticksize=0.05, col='red')							
    autocorr.plot(tmp[,1], auto.layout=FALSE, ask=FALSE, las=1, tcl=-0.2,
                  col='steelblue3', mgp=c(1,0.5,0))			
    
    traceplot(tmp[,4], col=c('blue','red'), lwd=0.5, lty=c(1,1), mgp=c(1,0.5,0), las=1,
              tcl=-0.2, main="", cex=0.7, smooth=TRUE)
    mtext(expression(Sigma[lat]), 2, 3)		
    foo1 = density(tmp[[1]][,4], bw = 1.06 * min(sd(tmp[[1]][,4]), IQR(tmp[[1]][,4])/1.34)
                   * length(tmp[[1]][,1])^-0.2)
    foo2 = density(tmp[[2]][,4], bw = 1.06 * min(sd(tmp[[2]][,4]), IQR(tmp[[2]][,4])/1.34)
                   * length(tmp[[2]][,4])^-0.2)			
    plot(foo1, col='blue', lwd=0.5, mgp=c(1,0.5,0), las=1, tcl=-0.2, main="",
         cex=0.7, ylim=c(-0.03*max(foo1$y,foo2$y), max(foo1$y,foo2$y)), 
         xlim=c(min(foo1$x,foo2$x),max(foo1$x,foo2$x)))
    lines(foo2, col='red', lwd=0.5)
    rug(tmp[[1]][,4], lwd=0.2, ticksize=0.05, col='blue')
    rug(tmp[[2]][,4], lwd=0.2, ticksize=0.05, col='red')				
    autocorr.plot(tmp[,4], auto.layout=FALSE, ask=FALSE, las=1, tcl=-0.2,
                  col='steelblue3', mgp=c(1,0.5,0))			
    
    #psi
    tmp = as.mcmc.list(m$mcmc$psi)
    traceplot(tmp, col=c('blue','red'), lwd=0.5, lty=c(1,1), mgp=c(1,0.5,0), las=1,
              tcl=-0.2, main="", cex=0.7, smooth=TRUE)
    mtext(expression(psi), 2, 3)			
    foo1 = density(tmp[[1]], bw = 1.06 * min(sd(tmp[[1]][,1]), IQR(tmp[[1]])/1.34) * 
                     length(tmp[[1]])^-0.2)
    foo2 = density(tmp[[2]], bw = 1.06 * min(sd(tmp[[2]][,1]), IQR(tmp[[2]])/1.34) * 
                     length(tmp[[2]])^-0.2)			
    plot(foo1, col='blue', lwd=0.5, mgp=c(1,0.5,0), las=1, tcl=-0.2, main="",
         cex=0.7, ylim=c(-0.03*max(foo1$y,foo2$y), max(foo1$y,foo2$y)), 
         xlim=c(min(foo1$x,foo2$x),max(foo1$x,foo2$x)))
    lines(foo2, col='red', lwd=0.5)	
    rug(tmp[[1]], lwd=0.2, ticksize=0.05, col='blue')
    rug(tmp[[2]], lwd=0.2, ticksize=0.05, col='red')						
    autocorr.plot(tmp, auto.layout=FALSE, ask=FALSE, las=1, tcl=-0.2,
                  col='steelblue3', mgp=c(1,0.5,0))
  }
  
lapply(fit, dostuff)	
}