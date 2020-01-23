#' Fits hierarchical state-space models to Argos data
#'
#' Takes output from \code{dat4jags}, sets up initial values, calls JAGS, and
#' aggregates results. Intended for internal use, called by \code{fit_ssm}.
#'
#' @param d structured data from \code{dat4jags} to be passed to JAGS
#' @param model the state-space model to be fit: hDCRW or hDCRWS
#' @param adapt number of samples in adaptation/burnin phase
#' @param samples number of posterior samples
#' @param thin thinning factor to reduce posterior sample autocorrelation
#' @param chains number of parallel McMC chains to run
#' @param span span
#' @return Returns a list of McMC samples from marginal posteriors and a
#' summary \code{data.frame} of mean and median position estimates.
#' @seealso Function to be called by \code{\link{fit_ssm}}.
#' @importFrom rjags jags.samples jags.model
#' @importFrom lubridate as_datetime
#' @importFrom msm rtnorm
#' @importFrom tibble data_frame as_data_frame
#' @export
hssm  <-
  function (d,
            model = "hDCRWS",
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
      
      with(
        dd,
        list(
          id = id,
          y = y,
          ts = ts,
          idx = idx,
          w = ws,
          itau2 = itau2,
          nu = nu,
          Nx = nrow(xs),
          Ny = nrow(dd$y),
          xs = xs,
          isigma2 = isigma2,
          rho = rho
        )
      )
    }
    prep <- lapply(d, ssm1)

    ## data
    N <- length(prep)
    Nx <- sapply(prep, function(x)
      x$Nx)
    Ny <- sapply(prep, function(x)
      x$Ny)
    y <- do.call(rbind, lapply(prep, function(x)
      x$y))
    itau2 <- do.call(rbind, lapply(prep, function(x)
      x$itau2))
    nu <- do.call(rbind, lapply(prep, function(x)
      x$nu))
    w <- as.numeric(unlist(lapply(prep, function(x)
      x$w)))
    idx <- lapply(prep, function(x)
      x$idx)
    m.idx <- cumsum(sapply(idx, function(x)
      max(x)))
    idx <-
      c(unlist(idx[[1]]), unlist(sapply(2:N, function(i)
        m.idx[i - 1] + idx[[i]])))
    
    Xidx <- cumsum(c(1, Nx))
    Yidx <- cumsum(c(1, Ny))
    data <-
      list(
        y = y,
        idx = idx,
        w = w,
        itau2 = itau2,
        nu = nu,
        Xidx = Xidx,
        Yidx = Yidx,
        N = N
      )
    
    ## inits
    isigma2 <- sapply(prep, function(x)
      x$isigma2)
    rho <- sapply(prep, function(x)
      x$rho)
    xs <- do.call(rbind, lapply(prep, function(x)
      x$xs))
    
    init.fn <- function(isigma2, rho, xs, N) {
      isigma2 <-
        rlnorm(2, log(c(mean(isigma2[1,]), mean(isigma2[2,]))), 0.1)
      rho <- rtnorm(1, mean(rho), 0.1, lower = -1, upper = 1)
      iSigma <- matrix(c(isigma2[1], rho, rho, isigma2[2]), 2, 2)
      gamma <- c(rbeta(1, 20, 20), NA)
      dev <- rbeta(1, 1, 1)
      alpha <- rbeta(2, 1, 1)
      lambda <- c(rbeta(1, 1, 1), NA)
      logpsi <- runif(N,-1, 1)
      x <-
        cbind(rnorm(nrow(xs), xs[, 1], 0.1), rnorm(nrow(xs), xs[, 2], 0.1))
      b <- rbinom(nrow(xs), 1, 0.5) + 1
      deviance <- rnorm(1)
      pD <- runif(1, 1000)
      
      init <-
        list(
          iSigma = iSigma,
          gamma = gamma[1],
          logpsi = logpsi,
          x = x,
          deviance = deviance,
          pD = pD
        )
      if (model == "hDCRWS") {
        init <-
          list(
            iSigma = iSigma,
            gamma = gamma,
            dev = dev,
            alpha = alpha,
            lambda = lambda,
            logpsi = logpsi,
            x = x,
            b = b,
            deviance = deviance,
            pD = pD
          )
      }
      init
    }
    inits <-
      lapply(1:chains, function(i)
        init.fn(isigma2, rho, xs, N))
    
    ## params
    params <- c("Sigma", "x", "gamma", "psi", "deviance", "pD")
    if (model == "hDCRWS") { params <- c(params, "alpha", "b") }

    model.file <-
      paste(system.file('jags', package = 'bsam'), "/", model, ".txt", sep = "")
    
    ## load DIC module in JAGS, if not already loaded
    ## this ensures deviance chains are monitored
    mods <- list.modules()
    if(!"dic" %in% mods) {
      load.module("dic")
    }
    
    burn <-
      jags.model(model.file,
                        data,
                        inits,
                        n.chains = chains,
                        n.adapt = adapt / 2
                 )
    update(burn, n.iter = adapt / 2)
    psamples <-
      jags.samples(burn, params, n.iter = samples, thin = thin)
    
    lon = apply(psamples$x[, 1, ,], 1, mean)
    lat = apply(psamples$x[, 2, ,], 1, mean)
    lon.q = apply(psamples$x[, 1, ,], 1, quantile, c(0.025, 0.5, 0.975))
    lat.q = apply(psamples$x[, 2, ,], 1, quantile, c(0.025, 0.5, 0.975))
   
    dts <-
      unlist(lapply(prep, function(x)
        x$ts))
    id  <- rep(as.character(sapply(prep, function(x)
      x$id)), Nx)
   
    summary <- data_frame(
      id = id,
      date = as_datetime(dts, tz = "GMT"),
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
    data <- do.call(rbind, lapply(d, function(x)
      x$obs))
    if (model == "hDCRWS") {
      b <- apply(psamples$b, 1, mean)
      b.5 <- apply(psamples$b, 1, median)
      summary <- as_data_frame(cbind(summary, b = b, b.5 = b.5))
    }
    list(
      summary = summary,
      mcmc = psamples,
      model = model,
      mcmc.settings = mcmc.settings,
      timestep = d$tstep,
      N = N,
      Nx = Nx,
      data = data,
      inits = inits
    )
  }
