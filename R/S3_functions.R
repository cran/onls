convInfo <- function(x, digits, show. = getOption("show.nls.convergence", TRUE)) 
{
  if (!is.null(x$convInfo)) 
    with(x$convInfo, {
      if (identical(x$call$algorithm, "port")) 
        cat("\nAlgorithm \"port\", convergence message: ", 
            stopMessage, "\n", sep = "")
      else {
        if (!isConv || show.) {
          cat("\nNumber of iterations", if (isConv) 
            "to convergence:"
            else "till stop:", finIter, "\nAchieved convergence tolerance:", 
            format(finTol, digits = digits))
          cat("\n")
        }
        if (!isConv) {
          cat("Reason stopped:", stopMessage)
          cat("\n")
        }
      }
    })
  invisible()
}

coef.onls <- function(object, ...)
{
  unlist(object$parONLS)
}

formula.onls <- function(x, ...)
{
  x$formula
}

weights.onls <- function(object, ...)
{
  object$weights
}

vcov.onls <- function(object, ...)
{
  ## object$vcov is already correctly scaled by onls() itself (see
  ## summary.onls()'s cov.unscaled logic / Details)
  object$vcov
}

df.residual.onls <- function(object, ...) object$df_resid 

residuals.onls <- function(object, ...)
{
  val <- as.vector(object$residONLS)
  if (!is.null(object$na.action)) 
    val <- napredict(object$na.action, val)
  lab <- "Vertical residuals from orthogonal model"
  attr(val, "label") <- lab
  val  
}

fitted.onls <- function(object, ...)
{
  val <- as.vector(object$fittedONLS)
  if (!is.null(object$na.action)) 
    val <- napredict(object$na.action, val)
  lab <- "Fitted values from orthogonal model"
  attr(val, "label") <- lab
  val  
}

deviance.onls <- function(object, ...)
{
  val <- sum(object$weights * object$residONLS^2)
  lab <- "(Weighted) Deviance (RSS) of vertical residuals from orthogonal model"
  attr(val, "label") <- lab
  val  
}

logLik.onls <- function(object, REML = FALSE, ...) 
{
  if (REML) stop("cannot calculate REML log-likelihood for \"onls\" objects")
  res <- object$residONLS
  N <- length(res)
  if (is.null(w <- object$weights)) w <- rep_len(1, N)
  
  ## Zero-weight observations are excluded entirely (matching R's own
  ## stats:::logLik.lm), rather than neutralized in place
  excl <- w == 0
  if (any(excl)) {
    res <- res[!excl]
    w   <- w[!excl]
    N   <- length(res)
  }
  
  val <- 0.5 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + log(sum(w * res^2))))
  
  q_free <- object$q - if (!is.null(object$fixed)) sum(as.logical(object$fixed)) else 0L
  attr(val, "df") <- 1L + q_free
  attr(val, "nobs") <- attr(val, "nall") <- N
  attr(val, "label") <- "Log-likelihood using vertical residuals from orthogonal model"
  class(val) <- "logLik.onls"
  val
}

print.logLik.onls <- function(x, digits = getOption("digits"), ...) 
{
  cat("'log Lik.' ", paste(format(c(x), digits = digits), collapse = ", "), 
      " (df=", format(attr(x, "df")), ")\n", 
      "\"", attr(x, "label"), "\"", sep = "")
  invisible(x)
}

predict.onls <- function(object, newdata, se.fit = FALSE, scale = NULL, df = NULL,
                         interval = c("none", "confidence", "prediction"),
                         level = 0.95, ...)
{
  interval <- match.arg(interval)
  xnms  <- object$pred_names
  pnms  <- object$param_names
  theta <- object$parONLS
  
  if (missing(newdata)) {
    newdata <- as.data.frame(object$X_obs)  
  } else {
    missing_cols <- setdiff(xnms, names(newdata))
    if (length(missing_cols)) stop("`newdata` is missing predictor column(s): ", paste(missing_cols, collapse = ", "))
    if (!is.null(cl <- object$dataClasses)) {
      cl_pred <- cl[intersect(names(cl), names(newdata))]
      if (length(cl_pred)) .checkMFClasses(cl_pred, newdata)
    }
  }
  
  env_base <- list()
  for (nm in xnms) env_base[[nm]] <- newdata[[nm]]
  for (nm in pnms) env_base[[nm]] <- theta[[nm]]
  env_base <- as.data.frame(env_base)
  
  fit <- eval(object$formula[[3L]], envir = env_base)
  if (!se.fit && interval == "none") return(fit)
  
  ## Only the FREE (non-fixed) parameters have rows/columns in
  ## object$vcov (fixed parameters carry zero uncertainty by
  ## construction, SE = NA/0 throughout the package)
  FIXSEL    <- if (!is.null(object$fixed)) which(as.logical(object$fixed)) else integer(0)
  freeSEL   <- setdiff(seq_along(theta), FIXSEL)
  pnms_free <- pnms[freeSEL]
  
  eps <- .Machine$double.eps^(1/3) * pmax(abs(unlist(theta)), 1e-3)
  J   <- matrix(NA_real_, length(fit), length(freeSEL), dimnames = list(NULL, pnms_free))
  for (jj in seq_along(freeSEL)) {
    j <- freeSEL[jj]
    env_up <- env_base
    env_up[[pnms[j]]] <- theta[[j]] + eps[j]
    env_dn <- env_base
    env_dn[[pnms[j]]] <- theta[[j]] - eps[j]
    f_up <- eval(object$formula[[3L]], envir = env_up)
    f_dn <- eval(object$formula[[3L]], envir = env_dn)
    J[, jj] <- (f_up - f_dn) / (2 * eps[j])
  }
  
  se_fit <- sqrt(pmax(rowSums((J %*% object$vcov) * J), 0))
  df_use <- if (is.null(df)) object$df_resid else df
  
  sigma_vert <- if (!is.null(scale)) scale else sqrt(sum(object$weights * object$residONLS^2) / object$df_resid)
  
  if (interval != "none") {
    tval <- qt(1 - (1 - level) / 2, df_use)
    se_used <- if (interval == "confidence") se_fit else sqrt(se_fit^2 + sigma_vert^2)
    out <- cbind(fit = fit, lwr = fit - tval * se_used, upr = fit + tval * se_used)
  } else out <- fit
  
  if (se.fit) return(list(fit = out, se.fit = se_fit, df = df_use, residual.scale = sigma_vert))
  out
}


print.onls <- function(x, ...) 
{
  digits <- max(3L, getOption("digits") - 3L)
  cat("Nonlinear orthogonal regression model\n")
  cat("  model: ", deparse(formula(x)), "\n", sep = "")
  cat("   data: ", deparse(x$data), "\n", sep = "")
  print(coef(x), digits = digits, ...)
  cat(" ", if (!is.null(x$weights) && diff(range(x$weights))) "weighted ", "vertical residual sum-of-squares: ", format(deviance(x), digits = digits), "\n", sep = "")
  cat(" ", if (!is.null(x$weights) && diff(range(x$weights))) "weighted ", "orthogonal residual sum-of-squares: ", format(deviance_o(x), digits = digits), "\n", sep = "")
  if (sum(x$ortho == TRUE) == length(x$ortho)) cat(" PASSED: ", sum(x$ortho == TRUE), " out of ", length(x$ortho), " fitted points are orthogonal.\n", sep = "")
  else cat(" FAILED: Only ", sum(x$ortho == TRUE), " out of ", length(x$ortho), " fitted points are orthogonal.\n", sep = "")
  convInfo(x, digits = digits)
  invisible(x)
}

summary.onls <- function(object, correlation = FALSE, symbolic.cor = FALSE, ...)
{
  residONLS <- as.vector(object$residONLS)
  param  <- coef(object)
  pnames <- names(param)
  q      <- length(pnames)
  
  FIXSEL  <- if (!is.null(object$fixed)) which(as.logical(object$fixed)) else integer(0)
  freeSEL <- setdiff(seq_len(q), FIXSEL)
  q_free  <- length(freeSEL)
  
  rdf <- object$df_resid
  
  resvarONLS <- if (rdf <= 0) NaN else deviance(object)/rdf
  resvar_o   <- if (rdf <= 0) NaN else deviance_o(object)/rdf
  
  Vcov <- object$vcov          # q_free x q_free
  se   <- object$std_errors    # length q, 0 at fixed positions
  
  tval <- pval <- setNames(rep(NA_real_, q), pnames)
  tval[freeSEL] <- param[freeSEL] / se[freeSEL]
  pval[freeSEL] <- 2 * pt(abs(tval[freeSEL]), rdf, lower.tail = FALSE)
  se_out <- se
  se_out[FIXSEL] <- NA_real_
  
  param <- cbind(param, se_out, tval, pval)
  dimnames(param) <- list(pnames, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  sigmaONLS <- sqrt(resvarONLS)
  
  chisq <- object$reduced_chisq
  cov.unscaled <- if (isTRUE(object$known_sigma) || is.null(chisq) ||
                      !is.finite(chisq) || chisq <= 0) Vcov else Vcov / chisq
  
  out <- list(formula = formula(object), residONLS = residONLS, sigmaONLS = sigmaONLS,
              sigma_o = sqrt(resvar_o), df = c(q_free, rdf), cov.unscaled = cov.unscaled,
              call = object$call, convInfo = object$convInfo, control = object$control,
              na.action = object$na.action, coefficients = param, parameters = param,
              known_sigma = object$known_sigma, reduced_chisq = object$reduced_chisq)
  
  if (correlation && rdf > 0) {
    se_free <- se[freeSEL]
    out$correlation <- Vcov / outer(se_free, se_free)
    dimnames(out$correlation) <- list(pnames[freeSEL], pnames[freeSEL])
    out$symbolic.cor <- symbolic.cor
  }
  class(out) <- "summary.onls"
  out
}

print.summary.onls <- function(x, digits = max(3L, getOption("digits") - 3L), 
                               symbolic.cor = x$symbolic.cor, 
                               signif.stars = getOption("show.signif.stars"), ...) 
{
  cat("\nFormula: ", paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n", sep = "")
  df <- x$df
  rdf <- df[2L]
  cat("\nParameters:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, ...)
  cat("\nResidual standard error of vertical distances:", format(signif(x$sigmaONLS, digits)), "on", rdf, "degrees of freedom")
  cat("\nResidual standard error of orthogonal distances:", format(signif(x$sigma_o, digits)), "on", rdf, "degrees of freedom")
  cat("\n")
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1) {
      cat("\nCorrelation of Parameter Estimates:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      }
      else {
        correl <- format(round(correl, 2), nsmall = 2L, 
                         digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }
  convInfo(x, digits = digits)
  if (nzchar(mess <- naprint(x$na.action))) cat("  (", mess, ")\n", sep = "")
  cat("\n")
  invisible(x)
}

plot.onls <- function(x, panel = NULL, fitted.nls = TRUE, fitted.onls = TRUE, asp = TRUE,
                      segments = TRUE, npoints = 200, nmesh = 25, ...)
{
  object  <- x
  xnms    <- object$pred_names
  pnms    <- object$param_names
  theta_o <- object$parONLS
  theta_n <- object$parNLS
  FORMULA <- object$formula[[3L]]
  nameY   <- attr(object$resp, "name")
  
  X  <- object$pred  
  Y  <- object$resp
  X0 <- object$x0
  Y0 <- object$y0
  
  dots  <- list(...)
  xlim0 <- dots$xlim; ylim0 <- dots$ylim
  
  p <- if (is.null(dim(X))) 1L else ncol(X)
  
  ## -----------------------------------------------------------------
  ## p == 2: interactive rgl 3D surface. Fully separate from everything
  ## below -- no single-predictor plotting logic applies here.
  ## -----------------------------------------------------------------
  if (p == 2L) {
    if (!requireNamespace("rgl", quietly = TRUE))
      stop("Package 'rgl' is required for a 2-predictor 3D plot; install it with ",
           "install.packages(\"rgl\").")
    
    x1nm <- xnms[1L]; x2nm <- xnms[2L]
    x1seq <- seq(min(X[, 1L], na.rm = TRUE), max(X[, 1L], na.rm = TRUE), length.out = nmesh)
    x2seq <- seq(min(X[, 2L], na.rm = TRUE), max(X[, 2L], na.rm = TRUE), length.out = nmesh)
    
    coefs <- object$coefficients
    grid  <- expand.grid(x1seq, x2seq)
    names(grid) <- c(x1nm, x2nm)
    for (nm in pnms) grid[[nm]] <- coefs[[nm]]
    zvals <- eval(FORMULA, envir = grid)
    zmat  <- matrix(zvals, nrow = nmesh, ncol = nmesh)
    
    rgl::open3d()
    rgl::surface3d(x1seq, x2seq, zmat, color = "lightblue", alpha = 0.6)
    rgl::surface3d(x1seq, x2seq, zmat, color = "grey20", alpha = 1,
                   front = "lines", back = "lines", lwd = 1.5, lit = FALSE)
    rgl::points3d(X[, 1L],  X[, 2L],  Y,  color = "black",   size = 6)
    rgl::points3d(X0[, 1L], X0[, 2L], Y0, color = "darkred", size = 5)
    
    if (segments) {
      seg_x <- as.vector(rbind(X[, 1L], X0[, 1L]))
      seg_y <- as.vector(rbind(X[, 2L], X0[, 2L]))
      seg_z <- as.vector(rbind(Y, Y0))
      rgl::segments3d(seg_x, seg_y, seg_z, color = "darkred", lwd = 1.5)
    }
    
    ## asp is now actually wired up here (it previously wasn't checked
    ## at all in the p == 2 branch): pass asp = FALSE to skip forcing a
    ## literal 1:1:1 view.
    if (isTRUE(asp)) rgl::aspect3d("iso")
    rgl::axes3d()
    rgl::box3d()
    rgl::title3d(xlab = x1nm, ylab = x2nm, zlab = nameY)
    
    return(invisible(NULL))
  }
  
  ## -----------------------------------------------------------------
  ## p == 1, or p > 2: build the list of single-predictor "plot specs"
  ## to draw -- one spec for p == 1 or a single chosen panel, or one
  ## spec PER PREDICTOR for the default multivariate grid (panel = NULL).
  ## -----------------------------------------------------------------
  if (p == 1L) {
    if (!is.null(panel)) warning("`panel` is ignored for single-predictor models.")
    specs <- list(list(xv = X, xname = xnms, x0v = X0, other_fixed = list(),
                       force_asp1 = TRUE, main_title = NULL))
  } else {
    reps <- setNames(colMeans(X0), xnms)
    
    if (!is.null(panel)) {
      if (is.character(panel)) {
        k <- panel[1L]
        if (!(k %in% xnms)) stop("`panel` must be one of: ", paste(xnms, collapse = ", "))
      } else {
        pidx <- suppressWarnings(as.integer(panel[1L]))
        if (is.na(pidx) || pidx < 1L || pidx > p)
          stop("`panel` must be an integer between 1 and ", p,
               ", or one of: ", paste(xnms, collapse = ", "))
        k <- xnms[pidx]
      }
      specs <- list(list(xv = X[, k], xname = k, x0v = X0[, k],
                         other_fixed = as.list(reps[setdiff(xnms, k)]),
                         force_asp1 = TRUE,
                         main_title = paste0(k, " (others fixed at mean foot value)")))
    } else {
      specs <- lapply(xnms, function(k) {
        list(xv = X[, k], xname = k, x0v = X0[, k],
             other_fixed = as.list(reps[setdiff(xnms, k)]),
             force_asp1 = FALSE,
             main_title = paste0(k, " (others fixed at mean foot value; see check_o() for orthogonality)"))
      })
    }
  }
  
  ## Grid layout only when drawing more than one spec (the default
  ## multivariate case with no panel chosen).
  if (length(specs) > 1L) {
    nr <- ceiling(sqrt(length(specs))); nc <- ceiling(length(specs) / nr)
    op <- par(mfrow = c(nr, nc)); on.exit(par(op))
  }
  
  ## -----------------------------------------------------------------
  ## Single shared plotting block, executed once per spec. This is the
  ## logic that previously lived in the separate .plot_one() /
  ## .square_companion() / .eval_curve() nested helper functions, now
  ## inlined directly here rather than as separate closures.
  ## -----------------------------------------------------------------
  for (spec in specs) {
    xv          <- spec$xv
    xname       <- spec$xname
    x0v         <- spec$x0v
    other_fixed <- spec$other_fixed
    force_asp1  <- spec$force_asp1
    main_title  <- spec$main_title
    
    plot_args <- dots
    if (isTRUE(asp) && force_asp1 && is.null(plot_args$asp)) plot_args$asp <- 1
    
    ## Companion xlim/ylim range: whichever of xlim/ylim the caller did
    ## NOT supply is derived (via extendrange()) from the data that
    ## falls inside the one that WAS supplied.
    xlim_req <- xlim0; ylim_req <- ylim0
    if (!is.null(xlim_req) && is.null(ylim_req)) {
      xlim_req <- sort(xlim_req)
      sel <- which(xv >= xlim_req[1L] & xv <= xlim_req[2L])
      yv_rel <- if (length(sel)) Y[sel] else Y
      ylim_req <- grDevices::extendrange(yv_rel)
    } else if (!is.null(ylim_req) && is.null(xlim_req)) {
      ylim_req <- sort(ylim_req)
      sel <- which(Y >= ylim_req[1L] & Y <= ylim_req[2L])
      xv_rel <- if (length(sel)) xv[sel] else xv
      xlim_req <- grDevices::extendrange(xv_rel)
    }
    plot_args$xlim <- xlim_req; plot_args$ylim <- ylim_req
    
    if (is.null(plot_args$xlab)) plot_args$xlab <- xname
    if (is.null(plot_args$ylab)) plot_args$ylab <- nameY
    if (is.null(plot_args$type)) plot_args$type <- "p"
    if (is.null(plot_args$pch))  plot_args$pch  <- 16
    if (is.null(plot_args$main)) plot_args$main <- main_title
    
    plot_args$x <- xv
    plot_args$y <- Y
    
    sweepX <- seq(min(xv, na.rm = TRUE), max(xv, na.rm = TRUE), length.out = npoints)
    
    ## Curve evaluation (previously .eval_curve()), inlined for each of
    ## the ONLS and NLS curves in turn.
    onlsY <- NULL
    if (fitted.onls) {
      env <- as.list(other_fixed)
      env[[xname]] <- sweepX
      for (nm in pnms) env[[nm]] <- theta_o[[nm]]
      onlsY <- eval(FORMULA, envir = as.data.frame(env))
    }
    nlsY <- NULL
    if (fitted.nls) {
      env <- as.list(other_fixed)
      env[[xname]] <- sweepX
      for (nm in pnms) env[[nm]] <- theta_n[[nm]]
      nlsY <- eval(FORMULA, envir = as.data.frame(env))
    }
    
    do.call(plot, plot_args)
    if (fitted.onls) lines(sweepX, onlsY, col = "darkred", lwd = 1.5)
    if (fitted.nls)  lines(sweepX, nlsY,  col = "blue", lwd = 1.5)
    if (segments) {
      points(x0v, Y0, col = "darkred", pch = 16, cex = 0.8)
      segments(x0v, Y0, xv, Y, col = "darkred", lwd = 1.5)
    }
    
    legend_labels <- c("ONLS fit", "NLS fit")[c(fitted.onls, fitted.nls)]
    legend_cols   <- c("darkred", "blue")[c(fitted.onls, fitted.nls)]
    if (length(legend_labels)) legend("topleft", legend = legend_labels, col = legend_cols,
                                      lty = 1, lwd = 2, bty = "n")
  }
  
  invisible(NULL)
}

profile.onls <- function(fitted, which = 1L:npar, maxpts = 100,
                         alphamax = 0.01, delta.t = cutoff/5, ...)
{
  f.summary <- summary(fitted)
  std.err <- f.summary$coefficients[, "Std. Error"]
  nobs <- length(resid(fitted))
  origPars <- unlist(fitted$parONLS)
  npar <- length(origPars)
  
  orig_fixed <- if (!is.null(fitted$fixed)) as.logical(fitted$fixed) else rep(FALSE, npar)
  
  lower <- fitted$lower
  lower <- rep_len(if (!is.null(lower)) as.double(lower) else -Inf, npar)
  upper <- fitted$upper
  upper <- rep_len(if (!is.null(upper)) as.double(upper) else Inf, npar)
  
  if (is.character(which))  which <- match(which, names(origPars), 0)
  which <- which[which >= 1 & which <= npar]
  which <- which[!orig_fixed[which]]   
  cutoff <- sqrt(qf(1 - alphamax, 1L, fitted$df_resid))
  out <- vector("list", npar)
  S.hat <- deviance_o(fitted)
  s2.hat <- deviance_o(fitted) / fitted$df_resid
  
  for (par in which) {
    pars <- origPars
    sgn <- -1
    count <- 1
    varying <- rep.int(TRUE, npar)
    varying[par] <- FALSE
    FIXED <- (!varying) | orig_fixed
    tau <- double(2 * maxpts)
    par.vals <- array(0, c(2L * maxpts, npar), list(NULL, names(pars)))
    tau[1L] <- 0
    par.vals[1, ] <- pars
    base <- pars[par]
    profile.par.inc <- delta.t * std.err[par]
    pars[par] <- base - profile.par.inc
    pars[par] <- pmin(upper[par], pmax(lower[par], pars[par]))
    
    while (count <= maxpts) {
      if (is.na(pars[par]) || isTRUE(all.equal(pars, par.vals[1, ]))
          || pars[par] < lower[par] || pars[par] > upper[par]
          || abs(pars[par] - base)/std.err[par] > 10 * cutoff) break
      
      FIT1 <- update(fitted, start = pars, fixed = FIXED, trace = FALSE)
      dev1 <- deviance_o(FIT1)
      fstat <- (dev1 - S.hat)/s2.hat
      
      if (is.na(fstat) || fstat < 0) break
      newtau <- sgn * sqrt(fstat)
      
      if (abs(newtau - tau[count]) < 0.1) break
      count <- count + 1
      tau[count] <- newtau
      par.vals[count, ] <- pars <- coef(FIT1)[1L:npar]
      
      if (abs(tau[count]) > cutoff) break
      pars <- pars + ((pars - par.vals[count - 1, ]) * delta.t)/abs(tau[count] - tau[count - 1])
      pars[-par] <- pmin(upper[-par], pmax(lower[-par], pars[-par]))
    }
    
    ind <- seq_len(count)
    tau[ind] <- tau[rev(ind)]
    par.vals[ind, ] <- par.vals[rev(ind), ]
    sgn <- 1
    newmax <- count + maxpts
    pars <- par.vals[count, ]
    pars[par] <- base + profile.par.inc
    pars[par] <- pmin(upper[par], pmax(lower[par], pars[par]))
    
    while (count <= newmax) {
      if (is.na(pars[par]) || isTRUE(all.equal(pars, par.vals[1, ]))
          || pars[par] < lower[par] || pars[par] >
          upper[par] || abs(pars[par] - base)/std.err[par] >
          10 * cutoff) break
      
      FIT2 <- update(fitted, start = pars, fixed = FIXED, trace = FALSE)
      dev2 <- deviance_o(FIT2)
      fstat2 <- (dev2 - S.hat)/s2.hat
      
      if (is.na(fstat2) || fstat2 < 0) break
      newtau <- sgn * sqrt(fstat2)
      if (abs(newtau - tau[count]) < 0.1) break
      count <- count + 1
      tau[count] <- newtau
      par.vals[count, ] <- pars <- coef(FIT2)[1L:npar]
      if (abs(tau[count]) > cutoff) break
      pars <- pars + ((pars - par.vals[count - 1, ]) * delta.t)/abs(tau[count] - tau[count - 1])
      pars[-par] <- pmin(upper[-par], pmax(lower[-par], pars[-par]))
    }
    
    ind <- seq_len(count)
    out[[par]] <- structure(list(tau = tau[ind], par.vals = par.vals[ind, , drop = FALSE]),
                            class = "data.frame", row.names = as.character(ind),
                            parameters = list(par = par, std.err = std.err[par]))
    
  }
  
  names(out)[which] <- names(coef(fitted))[which]
  out <- out[which]
  attr(out, "original.fit") <- fitted
  attr(out, "summary") <- f.summary
  class(out) <- c("profile.onls", "profile")
  out
}

confint.onls <- function(object, parm = "all", level = 0.95, k = 200, ...) 
{
  if (parm != "all") stop("Only 'all' is allowed for 'parm'.")
  n_success <- 0L
  alpha <- 1 - level
    
  DAT0 <- eval(object$data)
  xnms <- object$pred_names
  boot_coefs <- matrix(NA_real_, nrow = k, ncol = length(coef(object)))
  coef0 <- coef(object)
  se0 <- object$std_errors
    
  while (n_success < k) {
    idx   <- sample(1:nrow(DAT0), replace = TRUE)
    BOOT <- tryCatch(update(object, data = DAT0[idx, ], start = object$start, trace = FALSE, ...), error = function(e) NULL)
      
    if (any(abs(coef(BOOT) - coef0) > 20 * se0, na.rm = TRUE)) next # safeguard with 20x SE on any coefficient
      
    if (!is.null(BOOT) && BOOT$convergence == 0 && all(BOOT$ortho)) {
        n_success <- n_success + 1L
        boot_coefs[n_success, ] <- coef(BOOT)
        if (n_success %% 10 == 0) cat(n_success) else cat(".")
        if (n_success %% 50 == 0) cat("\n")
        flush.console()
      } else next
  }
  cat("\n")
    
  OUT <- apply(boot_coefs, 2, quantile, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  OUT <- t(OUT)
  colnames(OUT) <- c(paste0(100 * alpha / 2, "%"), paste0(100 * (1 - alpha / 2), "%"))
  rownames(OUT) <- names(coef0)  
  return(OUT)
}