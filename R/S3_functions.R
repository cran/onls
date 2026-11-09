convInfo <- function(x, digits, show. = getOption("show.nls.convergence", TRUE)) 
{
  if (!is.null(x$convInfo)) 
    with(x$convInfo, {
      if (!isConv || show.) {
        cat("\nNumber of iterations", if (isConv) "to convergence:" else "till stop:", finIter, 
            "\nConvergence tolerance (ftol):", format(finTol, digits = digits))
        cat("\n")
      }
      if (!isConv) {
        cat("Reason stopped:", stopMessage)
        cat("\n")
      }
    })
  invisible()
}

.napred <- function(object, val, label) {
  if (!is.null(object$na.action)) val <- napredict(object$na.action, val)
  attr(val, "label") <- label
  val
}
.q_free <- function(object) length(object$coefficients) - (if (!is.null(object$fixed)) sum(as.logical(object$fixed)) else 0L)

coef.onls        <- function(object, ...) object$coefficients
formula.onls     <- function(x, ...) x$formula
weights.onls     <- function(object, ...) object$weights
vcov.onls        <- function(object, ...) object$vcov       # already scaled by onls() (see summary.onls)
df.residual.onls <- function(object, ...) object$df_resid 

residuals.onls <- function(object, ...)
  .napred(object, as.vector(object$y_obs - object$fittedONLS), "Vertical residuals from orthogonal model")

fitted.onls <- function(object, ...)
  .napred(object, as.vector(object$fittedONLS), "Fitted values from orthogonal model")

deviance.onls <- function(object, ...)
{
  val <- sum(object$weights * (object$y_obs - object$fittedONLS)^2)
  attr(val, "label") <- "(Weighted) Deviance (RSS) of vertical residuals from orthogonal model"
  val  
}

logLik.onls <- function(object, REML = FALSE, ...) 
{
  if (REML) stop("cannot calculate REML log-likelihood for \"onls\" objects")
  res <- object$y_obs - object$fittedONLS
  w   <- object$weights
  keep <- w > 0                     # zero-weight observations are excluded, as in stats:::logLik.lm
  res <- res[keep]; w <- w[keep]; N <- length(res)
  val <- 0.5 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) + log(sum(w * res^2))))
  attr(val, "df") <- 1L + .q_free(object)
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

## Predictions at new predictor values. Confidence intervals: uncertainty of f(x, theta) from the
## parameters. Prediction intervals: a new observation (y, x) measured with the errors of the
## model, Var = Var(f) + s^2 (sigma_y^2 + g' Sigma_x g), g = df/dx at x (effective variance),
## with s^2 = 1 for known variances and the reduced chi-square otherwise.
## g' Sigma_x g for new observations, sigma_x in any form of onls() (standard deviations, or a
## p x p covariance matrix); zeros are allowed (error-free predictor)
.gSg_new <- function(G, sigma_x)
{
  m <- nrow(G); p <- ncol(G)
  if (anyNA(sigma_x)) stop("`sigma_x` must not contain NA.")
  if (is.matrix(sigma_x) && all(dim(sigma_x) == c(p, p))) {                ## a p x p matrix is always a covariance
    if (!isSymmetric(unname(sigma_x))) stop("A p x p `sigma_x` is a covariance matrix and must be symmetric.")
    return(pmax(.rowSums((G %*% sigma_x) * G, m, p), 0))
  }
  if (is.data.frame(sigma_x)) sigma_x <- as.matrix(sigma_x)              ## per-row SDs (also for m = p)
  if (any(sigma_x < 0)) stop("`sigma_x` must be non-negative.")
  S <- if (is.matrix(sigma_x)) sigma_x
       else if (length(sigma_x) == 1L || length(sigma_x) == p) matrix(sigma_x, m, p, byrow = TRUE)
       else if (p == 1L && length(sigma_x) == m) matrix(sigma_x, m, 1L)
       else stop("`sigma_x` must be a scalar, a length-p vector, a p x p covariance matrix or an m x p matrix.")
  .rowSums((G * S)^2, m, p)
}

predict.onls <- function(object, newdata, se.fit = FALSE, scale = NULL, df = NULL,
                         interval = c("none", "confidence", "prediction"),
                         level = 0.95, sigma_x = NULL, sigma_y = NULL, ...)
{
  interval <- match.arg(interval)
  xnms <- .xnames(object); pnms <- .pnames(object); theta <- object$coefficients
  if (missing(newdata)) newdata <- as.data.frame(object$X_obs)
  else {
    missing_cols <- setdiff(xnms, names(newdata))
    if (length(missing_cols)) stop("`newdata` is missing predictor column(s): ", paste(missing_cols, collapse = ", "))
    if (!is.null(cl <- object$dataClasses)) {
      cl_pred <- cl[intersect(names(cl), names(newdata))]
      if (length(cl_pred)) .checkMFClasses(cl_pred, newdata)
    }
  }
  Xnew  <- as.matrix(as.data.frame(newdata)[xnms]); m <- nrow(Xnew)
  env_f <- environment(object$formula); if (is.null(env_f)) env_f <- globalenv()
  rhs   <- object$formula[[3L]]
  f_eval <- function(th, Xi, strict = FALSE) .eval_nl(rhs, th, pnms, Xi, xnms, env_f, strict)
  fit <- f_eval(theta, Xnew, strict = TRUE)
  if (!se.fit && interval == "none") return(fit)
  
  gr   <- .make_model_grad(rhs, pnms, xnms, f_eval, env_f)(theta, Xnew)
  free <- if (!is.null(object$fixed)) !as.logical(object$fixed) else rep(TRUE, length(theta))
  J    <- gr$Jth[, free, drop = FALSE]
  V    <- object$vcov[free, free, drop = FALSE]
  se_fit <- sqrt(pmax(.rowSums((J %*% V) * J, m, ncol(J)), 0))
  known  <- isTRUE(object$known_sigma)
  df_use <- if (!is.null(df)) df else if (known) Inf else object$df_resid
  s2_fit <- if (known) 1 else object$reduced_chisq                        ## scale in vcov(object)
  s2     <- if (!is.null(scale)) scale^2 else s2_fit
  if (!is.null(scale)) se_fit <- se_fit * sqrt(s2 / s2_fit)                ## as predict.lm(scale =)
  
  out <- fit; se_used <- se_fit
  if (interval == "prediction") {
    if (!is.null(sigma_y) && !(length(sigma_y) %in% c(1L, m))) stop("`sigma_y` must have length 1 or nrow(newdata).")
    sy <- if (!is.null(sigma_y)) rep_len(sigma_y, m) else if (length(object$sigma_y) == 1L) rep(object$sigma_y, m) else
      stop("The fit used observation-specific `sigma_y`: supply `sigma_y` for `newdata`.")
    if (is.null(sigma_y) && diff(range(object$weights[object$weights > 0])) > 0)
      warning("The fit used non-constant `weights`: the prediction interval is for a new observation of weight 1 ",
              "(standard deviation `sigma_y`); supply `sigma_y` for newdata to set its response error.")
    G <- matrix(gr$Gx, m, length(xnms))
    gSg <- if (!is.null(sigma_x)) .gSg_new(G, sigma_x)                 ## sigma_x = 0 allowed here
           else if (!isTRUE(object$per_obs_x)) .rowSums(.Lt_solve(object$Lx[rep(1L, m), , , drop = FALSE], G)^2, m, ncol(G))
           else stop("The fit used observation-specific `sigma_x`: supply `sigma_x` for `newdata`.")
    se_used <- sqrt(se_fit^2 + s2 * (sy^2 + gSg))
  }
  if (interval != "none") {
    qq  <- if (is.finite(df_use)) qt(1 - (1 - level) / 2, df_use) else qnorm(1 - (1 - level) / 2)
    out <- cbind(fit = fit, lwr = fit - qq * se_used, upr = fit + qq * se_used)
  }
  if (se.fit) return(list(fit = out, se.fit = se_fit, df = df_use, residual.scale = sqrt(s2)))
  out
}

print.onls <- function(x, ...) 
{
  digits <- max(3L, getOption("digits") - 3L)
  wtd <- if (diff(range(x$weights))) "weighted " else ""
  cat("Nonlinear orthogonal regression model\n")
  cat("  model: ", deparse(formula(x)), "\n", sep = "")
  cat("   data: ", deparse(x$data), "\n", sep = "")
  print(coef(x), digits = digits, ...)
  cat(" ", wtd, "vertical residual sum-of-squares: ", format(deviance(x), digits = digits), "\n", sep = "")
  cat(" ", wtd, "orthogonal residual sum-of-squares: ", format(deviance_o(x), digits = digits), "\n", sep = "")
  n_ok <- sum(x$ortho %in% TRUE)
  if (all(is.na(x$ortho))) cat(" Orthogonality of the fitted points could not be checked.\n")
  else if (n_ok == length(x$ortho)) cat(" PASSED: ", n_ok, " out of ", length(x$ortho), " fitted points are orthogonal.\n", sep = "")
  else cat(" FAILED: Only ", n_ok, " out of ", length(x$ortho), " fitted points are orthogonal.\n", sep = "")
  convInfo(x, digits = digits)
  invisible(x)
}

summary.onls <- function(object, correlation = FALSE, symbolic.cor = FALSE, ...)
{
  param  <- coef(object); pnames <- names(param); q <- length(pnames)
  freeSEL <- if (!is.null(object$fixed)) which(!as.logical(object$fixed)) else seq_len(q)
  rdf  <- object$df_resid
  known <- isTRUE(object$known_sigma)
  Vcov <- object$vcov
  se   <- object$std_errors
  ## known variances: normal distribution (z), otherwise t with df_resid degrees of freedom
  tval <- pval <- se_out <- setNames(rep(NA_real_, q), pnames)
  se_out[freeSEL] <- se[freeSEL]
  tval[freeSEL]   <- param[freeSEL] / se[freeSEL]
  pval[freeSEL]   <- 2 * (if (known) pnorm(abs(tval[freeSEL]), lower.tail = FALSE) else pt(abs(tval[freeSEL]), rdf, lower.tail = FALSE))
  param <- cbind(param, se_out, tval, pval)
  dimnames(param) <- list(pnames, c("Estimate", "Std. Error", if (known) c("z value", "Pr(>|z|)") else c("t value", "Pr(>|t|)")))
  
  chisq <- object$reduced_chisq
  cov.unscaled <- if (known || is.null(chisq) || !is.finite(chisq) || chisq <= 0) Vcov else Vcov / chisq
  out <- list(formula = formula(object), residONLS = residuals(object),
              sigmaONLS = if (rdf <= 0) NaN else sqrt(deviance(object) / rdf),
              sigma_o = if (rdf <= 0) NaN else sqrt(deviance_o(object) / rdf),
              df = c(length(freeSEL), rdf), cov.unscaled = cov.unscaled,
              call = object$call, convInfo = object$convInfo, control = object$control,
              na.action = object$na.action, coefficients = param, parameters = param,
              known_sigma = known, reduced_chisq = chisq)
  if (correlation && rdf > 0) {
    se_free <- se[freeSEL]
    out$correlation <- Vcov[freeSEL, freeSEL, drop = FALSE] / outer(se_free, se_free)
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
  rdf <- x$df[2L]
  cat("\nParameters", if (isTRUE(x$known_sigma)) " (variances known: z values)", ":\n", sep = "")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, ...)
  cat("\nResidual standard error of vertical distances:", format(signif(x$sigmaONLS, digits)), "on", rdf, "degrees of freedom")
  cat("\nResidual standard error of orthogonal distances:", format(signif(x$sigma_o, digits)), "on", rdf, "degrees of freedom")
  cat("\n")
  correl <- x$correlation
  if (!is.null(correl) && (p <- NCOL(correl)) > 1) {
    cat("\nCorrelation of Parameter Estimates:\n")
    if (is.logical(symbolic.cor) && symbolic.cor) print(symnum(correl, abbr.colnames = NULL))
    else {
      correl <- format(round(correl, 2), nsmall = 2L, digits = digits)
      correl[!lower.tri(correl)] <- ""
      print(correl[-1, -p, drop = FALSE], quote = FALSE)
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
  xnms    <- .xnames(object)
  pnms    <- .pnames(object)
  theta_o <- object$coefficients
  theta_n <- object$parNLS
  FORMULA <- object$formula[[3L]]
  env_f   <- environment(object$formula); if (is.null(env_f)) env_f <- globalenv()
  nameY   <- .yname(object)
  
  X  <- if (ncol(object$X_obs) == 1L) object$X_obs[, 1L] else object$X_obs
  Y  <- object$y_obs
  X0 <- if (ncol(object$xi) == 1L) object$xi[, 1L] else object$xi
  Y0 <- object$yi
  
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
    zvals <- eval(FORMULA, envir = grid, enclos = env_f)
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
      onlsY <- eval(FORMULA, envir = as.data.frame(env), enclos = env_f)
    }
    nlsY <- NULL
    if (fitted.nls) {
      env <- as.list(other_fixed)
      env[[xname]] <- sweepX
      for (nm in pnms) env[[nm]] <- theta_n[[nm]]
      nlsY <- eval(FORMULA, envir = as.data.frame(env), enclos = env_f)
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

## Analysis of variance for nested onls models on the orthogonal objective S
## (the criterion the models minimize), not on the vertical residual sum of
## squares that anova.nls() would use through deviance().
##
## Estimated variances (known_sigma = FALSE): extra-sum-of-squares F test,
##   F = ((S_a - S_b) / (df_a - df_b)) / (S_b / df_b),  b = the larger model of a pair.
## Known variances (known_sigma = TRUE): S is -2 log-likelihood up to a constant,
##   so S_a - S_b is a likelihood-ratio statistic, chi-square with df_a - df_b.

anova.onls <- function(object, ...)
{
  models <- c(list(object), list(...))
  if (length(models) < 2L)
    stop("anova() for 'onls' needs at least two nested models fitted to the same data.")
  if (!all(vapply(models, inherits, logical(1L), what = "onls")))
    stop("All models must be fits of class 'onls'.")

  ## ---- comparability: same observations, same precisions, same variance assumption
  ref <- models[[1L]]
  for (m in models[-1L]) {
    if (m$n != ref$n || !isTRUE(all.equal(unname(m$y_obs), unname(ref$y_obs))))
      stop("The models were not fitted to the same observations (response or number of observations differ).")
    if (!isTRUE(all.equal(unname(m$Q_yy), unname(ref$Q_yy))))
      stop("The models differ in `weights` or `sigma_y`: their orthogonal objectives are not comparable.")
    common <- intersect(colnames(m$X_obs), colnames(ref$X_obs))
    if (!isTRUE(all.equal(unname(m$X_obs[, common, drop = FALSE]), unname(ref$X_obs[, common, drop = FALSE]))))
      stop("The models were not fitted to the same predictor values.")
    same_x <- setequal(colnames(m$X_obs), colnames(ref$X_obs))
    if (same_x) {
      ok <- isTRUE(all.equal(unname(m$Lx[, common, common, drop = FALSE]), unname(ref$Lx[, common, common, drop = FALSE])))
    } else {
      ## different predictor sets: only uncorrelated predictor errors can be compared
      offd <- function(L) { p <- dim(L)[2L]; p > 1L && any(vapply(seq_len(p), function(a) any(abs(L[, a, -a]) > 0), logical(1L))) }
      if (offd(m$Lx) || offd(ref$Lx))
        stop("Models with different predictors and correlated predictor errors (p x p `sigma_x`) cannot be compared.")
      dg <- function(L, k) L[, k, k]
      ok <- all(vapply(common, function(k) isTRUE(all.equal(dg(m$Lx, k), dg(ref$Lx, k))), logical(1L)))
    }
    if (!ok) stop("The models differ in `sigma_x`: their orthogonal objectives are not comparable.")
    if (!identical(isTRUE(m$known_sigma), isTRUE(ref$known_sigma)))
      stop("The models differ in `known_sigma`: compare models fitted with the same variance assumption.")
  }
  if (any(vapply(models, function(m) m$convergence != 0L, logical(1L))))
    warning("Not all models converged: the comparison may be misleading.")

  known <- isTRUE(ref$known_sigma)
  S  <- vapply(models, function(m) m$objective, numeric(1L))
  df <- vapply(models, function(m) as.numeric(m$df_resid), numeric(1L))
  k  <- length(models)

  ## ---- consecutive comparisons, as anova.nls(): model i against model i - 1
  dDf <- c(NA, -diff(df))
  dS  <- c(NA, -diff(S))
  stat <- pval <- rep(NA_real_, k)
  for (i in seq_len(k)[-1L]) {
    if (dDf[i] == 0) next
    big <- if (df[i] < df[i - 1L]) i else i - 1L             # the model with more parameters
    num <- (S[i - 1L] - S[i]) / (df[i - 1L] - df[i])          # change per degree of freedom (sign-free)
    if (known) {
      stat[i] <- num * abs(dDf[i])
      pval[i] <- stats::pchisq(stat[i], abs(dDf[i]), lower.tail = FALSE)
    } else {
      stat[i] <- num / (S[big] / df[big])
      pval[i] <- stats::pf(stat[i], abs(dDf[i]), df[big], lower.tail = FALSE)
    }
  }
  if (any(stat < 0, na.rm = TRUE))
    warning("A model with more parameters has a larger objective: the models are not nested, ",
            "or a fit did not reach its minimum.")

  TAB <- data.frame(df, S, dDf, dS, stat, pval)
  names(TAB) <- c("Res.Df", "Res.Obj", "Df", "Obj", if (known) "Chisq" else "F value",
                  if (known) "Pr(>Chi)" else "Pr(>F)")
  rownames(TAB) <- as.character(seq_len(k))
  forms <- vapply(models, function(m) paste(deparse(formula(m), width.cutoff = 500L), collapse = " "), character(1L))
  head <- c(paste0("Analysis of the orthogonal objective (", if (known) "known variances: likelihood-ratio chi-square test"
                   else "estimated variances: F test", ")\n"),
            paste0("Model ", format(seq_len(k)), ": ", forms, collapse = "\n"))
  structure(TAB, heading = head, class = c("anova", "data.frame"))
}
