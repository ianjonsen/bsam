#' Fits state-space models to Argos data
#'
#' Takes output from \code{dat4jags}, sets up initial values, calls JAGS, and
#' aggregates results. Intended for internal use, called by \code{fit_ssm}.
#'
#' @param d structured data from \code{dat4jags} to be passed to JAGS
#' @param model the state-space model to be fit: DCRW or DCRWS
#' @param adapt number of samples in adaptation/burnin phase
#' @param samples number of posterior samples
#' @param thin thinning factor to reduce posterior sample autocorrelation
#' @param chains number of parallel McMC chains to run
#' @param span span
#' @return Returns a list of McMC samples from marginal posteriors and a
#' summary \code{data.frame} of mean and median position estimates.
#' @seealso Function to be called by \code{\link{fit_ssm}}.
#' @importFrom rjags jags.samples jags.model
#' @importFrom msm rtnorm
#' @importFrom tibble data_frame
#' @export
ssm <-
  function (d,
            model = "DCRW",
            adapt,
            samples,
            thin,
            chains,
            span)
  {
    ssm1 <- function(dd) {
      gamma <- 0.5
      fit.lon <-
        loess(
          lon ~ as.numeric(date),
          data = dd$obs,
          span = span,
          na.action = "na.exclude",
          control = loess.control(surface = "direct")
        )
      fit.lat <-
        loess(
          lat ~ as.numeric(date),
          data = dd$obs,
          span = span,
          na.action = "na.exclude",
          control = loess.control(surface = "direct")
        )
      
      ## Predict track, increments and stochastic innovations
      xs <-
        cbind(predict(fit.lon, newdata = data.frame(date = as.numeric(dd$ts))),
              predict(fit.lat, newdata = data.frame(date = as.numeric(dd$ts))))
      
      if(dim(xs)[1] < 4) stop("\n\nCannot calculate initial values; \n  deployment length and time step imply too few states\n\n")
      ds <- xs[-1,] - xs[-nrow(xs),]
      es <- ds[-1,] - gamma * ds[-nrow(ds),]
      
      ## Estimate process variance components for initial values
      V <- cov(es)
      isigma2 <- diag(V) ^ -1
      rho <- V[1, 2] / prod(sqrt(isigma2))
      
      data <-
        with(dd,
             list(
               y = y,
               idx = idx,
               w = ws,
               itau2 = itau2,
               nu = nu,
               Nx = nrow(xs),
               Ny = nrow(dd$y)
             ))
      
      ## inits
      init.fn <- function() {
        isigma2 <- rlnorm(2, log(isigma2), 0.1)
        rho <- rtnorm(1, rho, 0.1, lower = -1, upper = 1)
        iSigma <- matrix(c(isigma2[1], rho, rho, isigma2[2]), 2, 2)
        gamma <- c(rbeta(1, 20, 20), NA)
        dev <- rbeta(1, 1, 1)
        alpha <- rbeta(2, 1, 1)
        lambda <- c(rbeta(1, 1, 1), NA)
        logpsi <- runif(1,-1, 1)
        x <-
          cbind(rnorm(nrow(xs), xs[, 1], 0.1), rnorm(nrow(xs), xs[, 2], 0.1))
        b <- rbinom(nrow(xs), 1, 0.5) + 1
        
        init <-
          list(
            iSigma = iSigma,
            gamma = gamma[1],
            logpsi = logpsi,
            x = x
          )
        if (model == "DCRWS") {
          init <-
            list(
              iSigma = iSigma,
              gamma = gamma,
              dev = dev,
              alpha = alpha,
              lambda = lambda,
              logpsi = logpsi,
              x = x,
              b = b
            )
        }
        init
      }
      inits <- lapply(1:chains, function(i)
        init.fn())
      params <- c("Sigma", "x", "gamma", "psi")
      if (model == "DCRWS")
        params <- c(params, "alpha", "b")
      
      model.file <-
        file.path(system.file("jags", package = "bsam"),
                  paste(model, ".txt", sep = ""))
      burn <-
        jags.model(model.file,
                   data,
                   inits,
                   n.chains = chains,
                   n.adapt = adapt / 2)
      update(burn, n.iter = adapt / 2)
      psamples <-
        jags.samples(burn, params, n.iter = samples, thin = thin)
      
      lon <- apply(psamples$x[, 1, ,], 1, mean)
      lat <- apply(psamples$x[, 2, ,], 1, mean)
      lon.q <-
        apply(psamples$x[, 1, ,], 1, quantile, c(0.025, 0.5, 0.975))
      lat.q <-
        apply(psamples$x[, 2, ,], 1, quantile, c(0.025, 0.5, 0.975))
      
      summary <- data_frame(
        id = dd$id,
        date = dd$ts,
        lon,
        lat,
        lon.025 = lon.q[1, ],
        lon.5 = lon.q[2, ],
        lon.975 = lon.q[3, ],
        lat.025 = lat.q[1, ],
        lat.5 = lat.q[2, ],
        lat.975 = lat.q[3, ]
      )
      model <- model
      mcmc.settings <-
        list(
          burnin = adapt,
          posterior.samples = samples,
          thinning = thin,
          n.chains = chains
        )
      data <- dd$obs
      
      if (model == "DCRWS") {
        b <- apply(psamples$b, 1, mean)
        b.5 <- apply(psamples$b, 1, median)
        summary <- data.frame(summary, b = b, b.5 = b.5)
      }
      
      out <-
        list(
          summary = summary,
          mcmc = psamples,
          model = model,
          mcmc.settings = mcmc.settings,
          timestep = dd$tstep,
          Nx = nrow(xs),
          data = data,
          inits = inits
        )
      
      out
    }
    
    lapply(d, ssm1)
  }
