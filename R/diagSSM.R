#' diagSSM
#' 
#' @param fit.in 
#'
#' @param save.to.pdf 
#' @importFrom coda autocorr.plot traceplot as.mcmc.list
#' @export
diagSSM = function(fit.in, save.to.pdf=FALSE){

if(names(fit.in)[1]=="summary") fit.in = list(fit.in)

'doplot' = function(fit){
	md = fit$model

	if(md=="DCRW" || md=="hDCRW"){	
		if(save.to.pdf){
			pdf(paste("diag", unique(fit$summary$id), "_params.pdf", sep=""), width=6.5, 
				height=10)
			}
		else{
			if(.Platform$OS.type=="windows") windows()
			else if(.Platform$OS.type=="unix"){
				quartz(w=6,h=6)	
				}	
			}		
		layout(matrix(1:20,5,4,byrow=TRUE), widths=c(3,3,3,3), heights=rep(1,5))
			par(mar=c(1,1,1,2), oma=c(3,4,2,0))
		
		tmp = as.mcmc.list(fit$mcmc$gamma)
		traceplot(tmp, col=c('blue','red'), lwd=0.5, lty=c(1,1), mgp=c(1,0.5,0), las=1,
			tcl=-0.2, cex=0.7, smooth=TRUE)
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
		mtext(expression(gamma), 2, 14)			
		coda::autocorr.plot(tmp, auto.layout=FALSE, ask=FALSE, las=1, tcl=-0.2,
			col='steelblue3', mgp=c(1,0.5,0))
		mtext("ACF: chain 1", 3, -1, outer=T, adj=0.635, cex=0.7)
		mtext("ACF: chain 2", 3, -1, outer=T, adj=0.91, cex=0.7)			
		mtext(paste(md, ": ", unique(fit$summary$id), sep=""), 3, 0.5, outer=TRUE)

		#theta
		tmp = as.mcmc.list(fit$mcmc$theta)
		traceplot(tmp, col=c('blue','red'), lwd=0.5, lty=c(1,1), mgp=c(1,0.5,0), las=1,
			tcl=-0.2, main="", cex=0.7, smooth=TRUE)
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
		mtext(expression(theta *"  ("*rad*")"), 2, 14)			
		autocorr.plot(tmp, auto.layout=FALSE, ask=FALSE, las=1, tcl=-0.2,
			col='steelblue3', mgp=c(1,0.5,0))			

		#Sigma
		tmp = as.mcmc.list(fit$mcmc$Sigma)
		traceplot(tmp[,1], col=c('blue','red'), lwd=0.5, lty=c(1,1), mgp=c(1,0.5,0), las=1,
			tcl=-0.2, main="", cex=0.7, smooth=TRUE)
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
		mtext(expression(Sigma[lon]), 2, 14)			
		autocorr.plot(tmp[,1], auto.layout=FALSE, ask=FALSE, las=1, tcl=-0.2,
			col='steelblue3', mgp=c(1,0.5,0))			

		traceplot(tmp[,4], col=c('blue','red'), lwd=0.5, lty=c(1,1), mgp=c(1,0.5,0), las=1,
			tcl=-0.2, main="", cex=0.7, smooth=TRUE)
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
		mtext(expression(Sigma[lat]), 2, 14)			
		autocorr.plot(tmp[,4], auto.layout=FALSE, ask=FALSE, las=1, tcl=-0.2,
			col='steelblue3', mgp=c(1,0.5,0))			

		# psi
		tmp = as.mcmc.list(fit$mcmc$psi)
		traceplot(tmp, col=c('blue','red'), lwd=0.5, lty=c(1,1), mgp=c(1,0.5,0), las=1,
			tcl=-0.2, main="", cex=0.7, smooth=TRUE)
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
		mtext(expression(psi), 2, 14)			
		autocorr.plot(tmp, auto.layout=FALSE, ask=FALSE, las=1, tcl=-0.2,
			col='steelblue3', mgp=c(1,0.5,0))
		if(save.to.pdf) dev.off()									
		}
		
	else if(md!="DCRW"){	
		if(save.to.pdf){
			if(md=="DCRWS"){
				pdf(paste("diag", unique(fit$summary$id), "_params.pdf", sep=""), width=6.5, 
					height=10)
				}
			else if(md=="hDCRWS"){
				pdf(paste("diag_", md, "_params.pdf", sep=""), width=6.5, 
					height=10)
				}
			}		
		else{
			if(.Platform$OS.type=="windows") windows()
			else if(.Platform$OS.type=="unix"){
				quartz(w=6,h=8)	
				}	
			}	
		if(md=="DCRWS"){
			layout(matrix(1:36,9,4,byrow=TRUE), widths=c(3,3,3,3), heights=rep(1,9))
			}
		else if(md=="hDCRW" || md=="hDCRWS"){
			layout(matrix(1:32,8,4,byrow=TRUE), widths=c(3,3,3,3), heights=rep(1,8))		
			}	
		par(mar=c(1,1,1,2), oma=c(3,4,2,0))			
		#gamma
		tmp = as.mcmc.list(fit$mcmc$gamma)
		traceplot(tmp[,1], col=c('blue','red'), lwd=0.5, lty=c(1,1), mgp=c(1,0.5,0), las=1,
			tcl=-0.2, cex=0.7, smooth=TRUE)
		mtext("Trace", 3, 0, cex=0.7)	
		foo1 = density(tmp[[1]][,1], bw = 1.06 * min(sd(tmp[[1]][,1]), IQR(tmp[[1]][,1])/1.34) 
			* length(tmp[[1]][,1])^-0.2)
		foo2 = density(tmp[[2]][,1], bw = 1.06 * min(sd(tmp[[2]][,1]), IQR(tmp[[2]][,1])/1.34)
			* length(tmp[[2]][,1])^-0.2)			
		plot(foo1, col='blue', lwd=0.5, mgp=c(1,0.5,0), las=1, tcl=-0.2,
			cex=0.7, ylim=c(-0.03*max(foo1$y,foo2$y), max(foo1$y,foo2$y)), 
			xlim=c(min(foo1$x,foo2$x),max(foo1$x,foo2$x)), main="")
		mtext("Density", 3, 0, cex=0.7)	
		lines(foo2, col='red', lwd=0.5)
		rug(tmp[[1]][,1], lwd=0.2, ticksize=0.05, col='blue')
		rug(tmp[[2]][,1], lwd=0.2, ticksize=0.05, col='red')							
		mtext(expression(gamma[1]), 2, 14)			
		autocorr.plot(tmp[,1], auto.layout=FALSE, ask=FALSE, las=1, tcl=-0.2,
			col='steelblue3', mgp=c(1,0.5,0))
		mtext("ACF: chain 1", 3, -1, outer=T, adj=0.635, cex=0.7)
		mtext("ACF: chain 2", 3, -1, outer=T, adj=0.91, cex=0.7)			
		if(md=="DCRWS"){
			mtext(paste(md, ": ", unique(fit$summary$id), sep=""), 3, 0.5, outer=TRUE)
			}
		else if(md=="hDCRW" || md=="hDCRWS"){
			ids = as.character(unique(fit$summary$id))			
			mtext(paste(md, ": ", paste(ids, sep="", collapse=", "),sep=""), 3, 0.5,
				outer=TRUE)
			}	
		traceplot(tmp[,2], col=c('blue','red'), lwd=0.5, lty=c(1,1), mgp=c(1,0.5,0), las=1,
			tcl=-0.2, main="", cex=0.7, smooth=TRUE)
		foo1 = density(tmp[[1]][,2], bw = 1.06 * min(sd(tmp[[1]][,2]), IQR(tmp[[1]][,2])/1.34) 
			* length(tmp[[1]][,2])^-0.2)
		foo2 = density(tmp[[2]][,2], bw = 1.06 * min(sd(tmp[[2]][,2]), IQR(tmp[[2]][,2])/1.34)
			* length(tmp[[2]][,2])^-0.2)			
		plot(foo1, col='blue', lwd=0.5, mgp=c(1,0.5,0), las=1, tcl=-0.2, main="",
			cex=0.7, ylim=c(-0.03*max(foo1$y,foo2$y), max(foo1$y,foo2$y)), 
			xlim=c(min(foo1$x,foo2$x),max(foo1$x,foo2$x)))
		lines(foo2, col='red', lwd=0.5)
		rug(tmp[[1]][,2], lwd=0.2, ticksize=0.05, col='blue')
		rug(tmp[[2]][,2], lwd=0.2, ticksize=0.05, col='red')						
		mtext(expression(gamma[2]), 2, 14)			
		autocorr.plot(tmp[,2], auto.layout=FALSE, ask=FALSE, las=1, tcl=-0.2,
			col='steelblue3', mgp=c(1,0.5,0))

		#theta
		tmp = as.mcmc.list(fit$mcmc$theta)	
		traceplot(tmp[,1], col=c('blue','red'), lwd=0.5, lty=c(1,1), mgp=c(1,0.5,0), las=1,
			tcl=-0.2, main="", cex=0.7, smooth=TRUE)
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
		mtext(expression(theta[1] *"  ("*rad*")"), 2, 14)			
		autocorr.plot(tmp[,1], auto.layout=FALSE, ask=FALSE, las=1, tcl=-0.2,
			col='steelblue3', mgp=c(1,0.5,0))

		traceplot(tmp[,2], col=c('blue','red'), lwd=0.5, lty=c(1,1), mgp=c(1,0.5,0), las=1,
			tcl=-0.2, main="", cex=0.7, smooth=TRUE)
		foo1 = density(tmp[[1]][,2], bw = 1.06 * min(sd(tmp[[1]][,2]), IQR(tmp[[1]][,2])/1.34) 
			* length(tmp[[1]][,2])^-0.2)
		foo2 = density(tmp[[2]][,2], bw = 1.06 * min(sd(tmp[[2]][,2]), IQR(tmp[[2]][,2])/1.34)
			* length(tmp[[2]][,2])^-0.2)			
		plot(foo1, col='blue', lwd=0.5, mgp=c(1,0.5,0), las=1, tcl=-0.2, main="",
			cex=0.7, ylim=c(-0.03*max(foo1$y,foo2$y), max(foo1$y,foo2$y)), 
			xlim=c(min(foo1$x,foo2$x),max(foo1$x,foo2$x)))
		lines(foo2, col='red', lwd=0.5)	
		rug(tmp[[1]][,2], lwd=0.2, ticksize=0.05, col='blue')
		rug(tmp[[2]][,2], lwd=0.2, ticksize=0.05, col='red')				
		mtext(expression(theta[2] *"  ("*rad*")"), 2, 14)			
		autocorr.plot(tmp[,2], auto.layout=FALSE, ask=FALSE, las=1, tcl=-0.2,
			col='steelblue3', mgp=c(1,0.5,0))

		#alpha
		tmp = as.mcmc.list(fit$mcmc$alpha)	
		traceplot(tmp[,1], col=c('blue','red'), lwd=0.5, lty=c(1,1), mgp=c(1,0.5,0), las=1,
			tcl=-0.2, main="", cex=0.7, smooth=TRUE)
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
		mtext(expression(alpha[1]), 2, 14)			
		autocorr.plot(tmp[,1], auto.layout=FALSE, ask=FALSE, las=1, tcl=-0.2,
			col='steelblue3', mgp=c(1,0.5,0))

		traceplot(tmp[,2], col=c('blue','red'), lwd=0.5, lty=c(1,1), mgp=c(1,0.5,0), las=1,
			tcl=-0.2, main="", cex=0.7, smooth=TRUE)
		foo1 = density(tmp[[1]][,2], bw = 1.06 * min(sd(tmp[[1]][,2]), IQR(tmp[[1]][,2])/1.34) 
			* length(tmp[[1]][,2])^-0.2)
		foo2 = density(tmp[[2]][,2], bw = 1.06 * min(sd(tmp[[2]][,2]), IQR(tmp[[2]][,2])/1.34)
			* length(tmp[[2]][,2])^-0.2)			
		plot(foo1, col='blue', lwd=0.5, mgp=c(1,0.5,0), las=1, tcl=-0.2, main="",
			cex=0.7, ylim=c(-0.03*max(foo1$y,foo2$y), max(foo1$y,foo2$y)), 
			xlim=c(min(foo1$x,foo2$x),max(foo1$x,foo2$x)))
		lines(foo2, col='red', lwd=0.5)	
		rug(tmp[[1]][,2], lwd=0.2, ticksize=0.05, col='blue')
		rug(tmp[[2]][,2], lwd=0.2, ticksize=0.05, col='red')						
		mtext(expression(alpha[2]), 2, 14)			
		autocorr.plot(tmp[,2], auto.layout=FALSE, ask=FALSE, las=1, tcl=-0.2,
			col='steelblue3', mgp=c(1,0.5,0))
			
		#Sigma
		tmp = as.mcmc.list(fit$mcmc$Sigma)
		traceplot(tmp[,1], col=c('blue','red'), lwd=0.5, lty=c(1,1), mgp=c(1,0.5,0), las=1,
			tcl=-0.2, main="", cex=0.7, smooth=TRUE)
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
		mtext(expression(Sigma[lon]), 2, 14)			
		autocorr.plot(tmp[,1], auto.layout=FALSE, ask=FALSE, las=1, tcl=-0.2,
			col='steelblue3', mgp=c(1,0.5,0))			

		traceplot(tmp[,4], col=c('blue','red'), lwd=0.5, lty=c(1,1), mgp=c(1,0.5,0), las=1,
			tcl=-0.2, main="", cex=0.7, smooth=TRUE)
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
		mtext(expression(Sigma[lat]), 2, 14)			
		autocorr.plot(tmp[,4], auto.layout=FALSE, ask=FALSE, las=1, tcl=-0.2,
			col='steelblue3', mgp=c(1,0.5,0))			
		if(md=="DCRWS"){
			# psi
			tmp = as.mcmc.list(fit$mcmc$psi)
			traceplot(tmp, col=c('blue','red'), lwd=0.5, lty=c(1,1), mgp=c(1,0.5,0), las=1,
				tcl=-0.2, main="", cex=0.7, smooth=TRUE)
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
			mtext(expression(psi), 2, 14)			
			autocorr.plot(tmp, auto.layout=FALSE, ask=FALSE, las=1, tcl=-0.2,
				col='steelblue3', mgp=c(1,0.5,0))
			if(save.to.pdf) dev.off()											
			}	

		else if(md=="hDCRW" || md=="hDCRWS"){
			if(save.to.pdf){ 
				dev.off()		
				pdf(paste("diag_", md, "_psi.pdf", sep=""), width=6, height=6)
				}
			else{
				if(.Platform$OS.type=="windows") windows()
				else if(.Platform$OS.type=="unix"){
					quartz(w=6,h=6)	
					}	
				}	
			# psi's
			tmp = as.mcmc.list(fit$mcmc$psi)
			dim.psi = ncol(tmp[[1]])
			layout(matrix(1:(4*dim.psi),dim.psi,4,byrow=TRUE), widths=c(3,3,3,3),
				heights=rep(1,dim.psi))			
			par(mar=c(1,1,1,2), oma=c(3,4,2,0))			
			
			for(i in 1:dim.psi){
				traceplot(tmp[,i], col=c('blue','red'), lwd=0.5, lty=c(1,1), mgp=c(1,0.5,0),
					las=1, tcl=-0.2, cex=0.7, smooth=TRUE)
				if(i==1) mtext("Trace", 3, 0, cex=0.7)	
				foo1 = density(tmp[[1]][,i], bw = 1.06 * min(sd(tmp[[1]][,i]), 
					IQR(tmp[[1]][,i])/1.34) * length(tmp[[1]][,i])^-0.2)
				foo2 = density(tmp[[2]][,i], bw = 1.06 * min(sd(tmp[[2]][,i]), 
					IQR(tmp[[2]][,i])/1.34) * length(tmp[[2]][,i])^-0.2)			
				plot(foo1, col='blue', lwd=0.5, mgp=c(1,0.5,0), las=1, tcl=-0.2, main="",
					cex=0.7, ylim=c(-0.03*max(foo1$y,foo2$y), max(foo1$y,foo2$y)), 
					xlim=c(min(foo1$x,foo2$x),max(foo1$x,foo2$x)))
				if(i==1) mtext("Density", 3, 0, cex=0.7)	
				lines(foo2, col='red', lwd=0.5)	
				rug(tmp[[1]][,i], lwd=0.2, ticksize=0.05, col='blue')
				rug(tmp[[2]][,i], lwd=0.2, ticksize=0.05, col='red')						
				mtext(bquote(psi[.(i)]), 2, 14)	
				if(i==1){
					mtext("ACF: chain 1", 3, -1, outer=T, adj=0.635, cex=0.7)
					mtext("ACF: chain 2", 3, -1, outer=T, adj=0.91, cex=0.7)
					ids = as.character(unique(fit$summary$id))			
					mtext(paste(md, ": ", paste(ids, sep="", collapse=", "),sep=""), 3, 0.5,
						outer=TRUE)
					}			
				autocorr.plot(tmp[,i], auto.layout=FALSE, ask=FALSE, las=1, tcl=-0.2,
					col='steelblue3', mgp=c(1,0.5,0))			
				}
			if(save.to.pdf) dev.off() 
			}				
		}	

	}
	
lapply(fit.in, doplot)	
}