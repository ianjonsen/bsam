#' plotSSM
#' 
#' @param fit.in fit.in
#'
#' @param save.to.pdf save.to.pdf
#'
#' @export
#' @importFrom PBSmapping plotMap
`plotSSM` = function(fit.in, save.to.pdf=FALSE) 
{
options(warn=-1)
  on.exit(options(warn=0))
  worldLLhigh <- NULL
data("worldLLhigh", envir = environment())
cols = colorRampPalette(c("blue", "yellow"))
ncols = 21

'doplot' = function(fit){
	if(fit$model=="hDCRWS"){
		est = split(fit$summary, fit$summary$id)
		dat = split(fit$data, fit$data$id)
		Nest = sapply(est, nrow)
		idx = c(1,cumsum(Nest))
		mc = list()
		for(i in 1:length(Nest)){
			b = fit$mcmc$b[idx[i]:idx[i+1],,]	
			x = fit$mcmc$x[idx[i]:idx[i+1],,,]
			mc[[i]] = list(b=b, x=x)
			}
		par.means = list(Sigma=rbind(apply(fit$mcmc$Sigma[1,,,],1,mean), 
			apply(fit$mcmc$Sigma[2,,,],1,mean)), alpha = apply(fit$mcmc$alpha,1,mean),
			gamma = apply(fit$mcmc$gamma,1,mean), theta = apply(fit$mcmc$theta,1,mean))
		
		## do each plot for hDCRWS in it's own window
		sapply(1:length(Nest), function(i){
			xl1 = extendrange(dat[[i]]$lon, f=0.1) + 360
			xl2 = extendrange(dat[[i]]$lon, f=0.1)
			yl = extendrange(dat[[i]]$lat, f=0.1)
			if(save.to.pdf){
				pdf(paste("ssm", unique(est[[i]]$id), ".pdf", sep=""), width=6,
					height=5, pointsize=12)
				}
			else {
				if(.Platform$OS.type=="windows") windows()
				else if(.Platform$OS.type=="unix"){
					quartz(width = 6, height = 5)	
					}	
				}		
			layout(matrix(c(1,1,0),3,1,byrow=TRUE), widths=3, heights=c(4,4,1))	
			plotMap(worldLLhigh, xlim=xl1, ylim=yl, border=NA, bg=grey(0.7),
				col=grey(0.05), axes=FALSE, xlab="", ylab="")
			par(new=TRUE)	
	    	plotMap(worldLLhigh, xlim=xl2, ylim=yl, border=NA, col=grey(0.05), 
	    		axes=FALSE, xlab="", ylab="")    		
			axis(1, at=pretty(xl2, n=5), mgp=c(1,0.3,0),tcl=-0.3, cex.axis=0.8, lwd=0.5)
			axis(2, at=pretty(yl, n=5), mgp=c(1,0.3,0),tcl=-0.3, cex.axis=0.8, lwd=0.5,
				las=1)
			axis(3, at=pretty(xl2, n=5), labels = FALSE, mgp=c(1,0.3,0), tcl=-0.3, lwd=0.5)
			axis(4, at=pretty(yl, n=5), labels = FALSE, mgp=c(1,0.3,0), tcl=-0.3, lwd=0.5)						
		   	points(lat ~ lon, data = dat[[i]], type = "b", pch="+", col = grey(0.5))
			matpoints(mc[[i]]$x[,1,,], mc[[i]]$x[,2,,], 
				col = rgb(0.95, 0.95, 0.95, 0.3), lwd = 0.1, pch = ".")    	
			lines(lat ~ lon, est[[i]], col = rgb(0,0,1,0.75), lwd = 0.5)
	        points(lat ~ lon, data = est[[i]], pch = 21,
    	    	bg=cols(ncols)[floor((b-1)*ncols)+1], col = grey(0.9), lwd=0.5)   	    
    	    mtext(paste("Track id: ", unique(est[[i]]$id), "   Argos positions: ",
		        nrow(dat[[i]]), "   Estimated positions: ", nrow(est[[i]]), 
    		    "   Time step: ", with(est[[i]], difftime(date[2], date[1], 
    		    units = "days")), " days", sep=""), 1, -5, adj = 0.6, outer = TRUE,
    		 	cex=0.8)
    		mtext(paste("Mean turn angles (theta): ", 	
    	    	round(par.means$theta[1] * 180/pi, 2), ", ", 
    	    	round(par.means$theta[2] * 180/pi, 2),
    	    	"   Move persistence (gamma): ", round(par.means$gamma[1], 2), ", ",
    	    	round(par.means$gamma[2], 2), sep=""), 1, -3.5, adj = 0.5, outer = TRUE,
    	    	cex=0.8) 	
		    mtext(paste("Process variance (Sigma: lon): ", 	
		    	round(par.means$Sigma[1,1], 3), "   Process variance (Sigma: lat): ",
		    	round(par.means$Sigma[2,2], 3), sep=""), 1, -2, adj = 0.5, outer = TRUE,
		    	cex=0.8)  
		    mtext("Longitude", 1, 1.5, cex=0.9)
			mtext("Latitude", 2, 1.5, cex=0.9)  
			mtext(paste("Model: ", fit$model), 3, 1, adj=0.5, cex=0.9)	
			if(save.to.pdf) dev.off()    		
			})
		}	
	else if(fit$model!="hDCRWS"){		
		xl1 = extendrange(fit$data$lon, f=0.1) + 360
		xl2 = extendrange(fit$data$lon, f=0.1)
		yl = extendrange(fit$data$lat, f=0.1)
		if(save.to.pdf){
			pdf(paste("ssm", unique(fit$summary$id), ".pdf", sep=""), width=6,
				height=5, pointsize=12)
			}
		else {
			if(.Platform$OS.type=="windows") windows()
			else if(.Platform$OS.type=="unix"){
				quartz(width = 6, height = 5)	
				}
			}	
		layout(matrix(c(1,1,0),3,1,byrow=TRUE), widths=3, heights=c(4,4,1))	
		plotMap(worldLLhigh, xlim=xl1, ylim=yl, border=NA, bg=grey(0.7), col=grey(0.05),
		    axes=FALSE, xlab="", ylab="")
		par(new=TRUE)	
	    plotMap(worldLLhigh, xlim=xl2, ylim=yl, border=NA, col=grey(0.05), 
	    	axes=FALSE, xlab="", ylab="")    		
		axis(1, at=pretty(xl2, n=5), mgp=c(1,0.3,0),tcl=-0.3, cex.axis=0.8, lwd=0.5)
		axis(2, at=pretty(yl, n=5), mgp=c(1,0.3,0),tcl=-0.3, cex.axis=0.8, lwd=0.5,
			las=1)
		axis(3, at=pretty(xl2, n=5), labels = FALSE, mgp=c(1,0.3,0), tcl=-0.3, lwd=0.5)
		axis(4, at=pretty(yl, n=5), labels = FALSE, mgp=c(1,0.3,0), tcl=-0.3, lwd=0.5)						
	   	points(lat ~ lon, data = fit$data, type = "b", pch="+", col = grey(0.5))
		matpoints(fit$mcmc$x[,1,,], fit$mcmc$x[,2,,], 
			col = rgb(0.95, 0.95, 0.95, 0.3), lwd = 0.1, pch = ".")    	
		lines(lat ~ lon, fit$summary, col = rgb(0,0,1,0.75), lwd = 0.5)
		mtext(paste("Track id: ", unique(fit$summary$id), "   Argos positions: ", 
	        nrow(fit$data), "   Estimated positions: ", nrow(fit$summary), 
	        "   Time step: ", with(fit$summary, difftime(date[2], date[1], 
	        units = "days")), " days", sep=""), 1, -5, adj = 0.6, outer = TRUE, cex=0.8)
	    mtext("Longitude", 1, 1.5, cex=0.9)
		mtext("Latitude", 2, 1.5, cex=0.9)    
	    switch(fit$model, DCRW = {
		    points(lat~lon, data=fit$summary, pch=21, bg='blue', col=grey(0.9), lwd=0.5)
	        mtext(paste("Mean turn angle: ", 
	        	round(mean(fit$mcmc$theta) * 180/pi, 2), 
		        "   Move persistence: ",
		        round(mean(fit$mcmc$gamma), 2), sep=""), 1, -3.5, 
		        adj = 0.5, outer = TRUE, cex=0.8)
	    	}, DCRWS = {
        	points(lat ~ lon, data = fit$summary, pch = 21, 
        		bg=cols(ncols)[floor((b-1)*ncols)+1], col = grey(0.9), lwd=0.5)
    	    mtext(paste("Mean turn angles (theta): ", 	
    	    	round(mean(fit$mcmc$theta[1,,]) * 180/pi, 2), ", ", 
    	    	round(mean(fit$mcmc$theta[2,,]) * 180/pi, 2),
    	    	"   Move persistence (gamma): ", round(mean(fit$mcmc$gamma[1,,]), 2), 
    	    	", ", round(mean(fit$mcmc$gamma[2,,]), 2), sep=""), 1, -3.5, 
    	    	adj = 0.5, outer = TRUE, cex=0.8)
	    	})
	    mtext(paste("Process variance (Sigma: lon): ",
	    	round(mean(fit$mcmc$Sigma[1,1,,]),3), "   Process variance (Sigma: lat): ", 
	    	round(mean(fit$mcmc$Sigma[2,2,,]),3), sep=""), 1, -2, adj = 0.5, 
	    	outer = TRUE, cex=0.8)	
	    mtext(paste("Model: ", fit$model), 3, 1, adj=0.5, cex=0.9)		
	   	if(save.to.pdf) dev.off()	 	
		}
	}

if(names(fit.in)[1]=="summary"){
	input = list(fit.in)
	lapply(input, doplot)
	}	

else{
	lapply(fit.in, doplot)
	}

cat("\n Use diagSSM() to examine convergence criteria \n")
#options(warn=0)
}