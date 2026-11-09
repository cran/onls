## =============================================================================
## Refitting, profiles and confidence intervals for 'onls' objects
## =============================================================================

## Refit an 'onls' model from the data stored in the object (after `subset` and
## `na.action`), optionally on a row index `idx` (bootstrap) and with other
## `start`/`fixed` values (profiling). Row-aligned `weights`, `sigma_x` and
## `sigma_y` are taken with the same rows, so they stay paired with their data,
## and nothing has to be re-evaluated in the caller's environment.
.refit <- function(object, idx = NULL, start = coef(object), fixed = object$fixed) {
  n <- object$n
  i <- if (is.null(idx)) seq_len(n) else idx
  dat <- as.data.frame(object$X_obs[i, , drop = FALSE])
  dat[[.yname(object)]] <- object$y_obs[i]

  w  <- object$weights
  w  <- if (is.null(w) || all(w == 1)) NULL else w[i]
  sx <- object$sigma_x
  if (isTRUE(object$per_obs_x)) sx <- if (is.matrix(sx)) sx[i, , drop = FALSE] else sx[i]
  sy <- object$sigma_y
  if (length(sy) == n && n > 1L) sy <- sy[i]

  args <- list(formula = object$formula, data = dat, start = as.list(start),
               weights = w, sigma_x = sx, sigma_y = sy, known_sigma = object$known_sigma,
               extend = object$extend, window = object$window,
               control = if (is.null(object$control)) list() else object$control,
               lower = object$lower, upper = object$upper, fixed = fixed, trace = FALSE,
               algorithm = if (is.null(object$algorithm)) "auto" else object$algorithm)
  suppressWarnings(do.call(onls, args))
}

## Profile of the orthogonal objective S for each free parameter: the parameter is
## held fixed on a grid around its estimate, all other free parameters (and all foot
## points) are re-optimized, and tau = sign * sqrt((S - S_hat) / s2) is recorded, with
## s2 = S_hat / df_resid (estimated variances) or 1 (known variances).
profile.onls <- function(fitted, which = NULL, maxpts = 100, alphamax = 0.01,
                         delta.t = cutoff / 5, ...)
{
  cf  <- coef(fitted); pn <- names(cf); npar <- length(cf)
  fixed0 <- if (!is.null(fitted$fixed)) as.logical(fitted$fixed) else rep(FALSE, npar)
  se  <- fitted$std_errors
  lower <- rep_len(if (!is.null(fitted$lower)) as.double(fitted$lower) else -Inf, npar)
  upper <- rep_len(if (!is.null(fitted$upper)) as.double(fitted$upper) else  Inf, npar)

  if (is.null(which)) which <- which(!fixed0)
  if (is.character(which)) which <- match(which, pn, 0L)
  which <- which[which >= 1L & which <= npar]
  if (any(fixed0[which])) warning("Fixed parameters are not profiled: ", paste(pn[which][fixed0[which]], collapse = ", "))
  which <- which[!fixed0[which]]
  if (!length(which)) stop("No free parameters to profile.")

  known  <- isTRUE(fitted$known_sigma)
  S.hat  <- fitted$objective
  if (!known && (sum(fitted$Q_yy > 0) - sum(!fixed0) <= 0 || S.hat <= 1e-12 * sum(fitted$y_obs^2)))
    stop("No residual degrees of freedom (exact fit): profiles of the objective are not defined.")
  s2.hat <- if (known) 1 else S.hat / fitted$df_resid
  cutoff <- if (known) sqrt(qchisq(1 - alphamax, 1L)) else sqrt(qf(1 - alphamax, 1L, fitted$df_resid))
  better <- FALSE
  
  ## one direction (sgn = -1 / +1) of the profile of parameter j; `reason` says why it ended
  walk <- function(j, sgn, FIX) {
    tau <- 0; pv <- matrix(cf, 1L, npar, dimnames = list(NULL, pn)); reason <- "maxpts"
    pars <- cf; pars[j] <- cf[j] + sgn * delta.t * se[j]
    for (k in seq_len(maxpts)) {
      if (!is.finite(pars[j])) { reason <- "range"; break }
      at_bound <- FALSE
      if (pars[j] < lower[j] || pars[j] > upper[j]) {       # evaluate at the bound, then stop
        bnd <- if (pars[j] < lower[j]) lower[j] else upper[j]
        if (pv[nrow(pv), j] == bnd) { reason <- "bound"; break }
        pars[j] <- bnd; at_bound <- TRUE
      }
      if (abs(pars[j] - cf[j]) / se[j] > 10 * cutoff) { reason <- "range"; break }
      fit <- tryCatch(.refit(fitted, start = pars, fixed = FIX), error = function(e) e)
      if (inherits(fit, "error")) { reason <- paste("a refit failed:", conditionMessage(fit)); break }
      if (fit$convergence != 0L) { reason <- "a refit did not converge"; break }
      fstat <- (fit$objective - S.hat) / s2.hat
      if (!is.finite(fstat)) { reason <- "a refit gave a non-finite objective"; break }
      if (fstat < -1e-8 * max(1, S.hat / s2.hat)) { better <<- TRUE; reason <- "better"; break }
      newtau <- sgn * sqrt(max(fstat, 0))
      if (abs(newtau - tau[length(tau)]) < 0.1 && !at_bound) { reason <- "flat"; break }
      prev <- pv[nrow(pv), ]
      tau <- c(tau, newtau); pv <- rbind(pv, coef(fit))
      if (abs(newtau) > cutoff) { reason <- "cutoff"; break }
      if (at_bound) { reason <- "bound"; break }
      pars <- coef(fit) + (coef(fit) - prev) * delta.t / abs(newtau - tau[length(tau) - 1L])
      pars[-j] <- pmin(upper[-j], pmax(lower[-j], pars[-j]))
    }
    list(tau = tau[-1L], pv = pv[-1L, , drop = FALSE], reason = reason)
  }
  
  out <- vector("list", length(which)); names(out) <- pn[which]
  for (j in which) {
    if (!is.finite(se[j]) || se[j] <= 0) { warning("No standard error for ", pn[j], ": not profiled."); next }
    FIX <- fixed0; FIX[j] <- TRUE
    dn <- walk(j, -1, FIX); up <- walk(j, 1, FIX)
    ord <- rev(seq_along(dn$tau))
    tau <- c(dn$tau[ord], 0, up$tau)
    pv  <- rbind(dn$pv[ord, , drop = FALSE], cf, up$pv)
    rownames(pv) <- NULL
    out[[pn[j]]] <- structure(list(tau = tau, par.vals = pv), class = "data.frame",
                              row.names = as.character(seq_along(tau)),
                              parameters = list(par = j, std.err = se[j],
                                                stop = c(lower = dn$reason, upper = up$reason)))
  }
  out <- out[!vapply(out, is.null, logical(1L))]
  if (better) warning("Profiling found a lower objective than the fit: the original fit is not at the minimum.")
  attr(out, "original.fit") <- fitted
  attr(out, "summary") <- summary(fitted)
  class(out) <- c("profile.onls", "profile")
  out
}

## quantile used for intervals: normal for known variances, t otherwise
.ci_quantile <- function(object, level) {
  a <- (1 - level) / 2
  if (isTRUE(object$known_sigma)) qnorm(1 - a) else qt(1 - a, object$df_resid)
}

.ci_names <- function(level) {
  a <- (1 - level) / 2
  paste(format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE, digits = 3), "%")
}

.why_na <- function(reason) {
  switch(sub(":.*", "", reason),
         flat   = "the profile flattens out, so the data do not bound the parameter on this side",
         bound  = "the limit lies beyond the `lower`/`upper` bound of the parameter",
         cutoff = "the profile was computed for a lower level; use a smaller `alphamax`",
         maxpts = "the profile stopped after `maxpts` points; try a larger `maxpts`",
         range  = "the profile ran 10 cut-offs (in standard errors) from the estimate without reaching the level, so the data hardly bound the parameter on this side",
         better = "a refit found a lower objective than the fit (the fit is not at the minimum)",
         nonmonotone = "the profile is not monotone (refits reached another minimum)",
         reason)
}

confint.profile.onls <- function(object, parm, level = 0.95, ...)
{
  fit <- attr(object, "original.fit")
  if (missing(parm)) parm <- names(object)
  if (is.numeric(parm)) parm <- names(coef(fit))[parm]
  parm <- intersect(parm, names(object))
  cut  <- .ci_quantile(fit, level)
  out  <- matrix(NA_real_, length(parm), 2L, dimnames = list(parm, .ci_names(level)))
  msgs <- character(0)
  for (nm in parm) {
    d   <- object[[nm]]; tau <- d$tau; pv <- d$par.vals[, nm]
    why <- attr(d, "parameters")$stop
    i0  <- which(tau == 0)[1L]
    ## monotone stretch on each side of the estimate
    lo <- i0; while (lo > 1L && is.finite(tau[lo - 1L]) && tau[lo - 1L] < tau[lo] && pv[lo - 1L] < pv[lo]) lo <- lo - 1L
    hi <- i0; while (hi < length(tau) && is.finite(tau[hi + 1L]) && tau[hi + 1L] > tau[hi] && pv[hi + 1L] > pv[hi]) hi <- hi + 1L
    if (lo > 1L) why["lower"] <- "nonmonotone"
    if (hi < length(tau)) why["upper"] <- "nonmonotone"
    tt <- tau[lo:hi]; xx <- pv[lo:hi]
    if (length(tt) >= 2L) {
      inv <- tryCatch({                                  # monotone spline of tau over the parameter, inverted
        f  <- splinefun(xx, tt, method = "hyman")
        xs <- seq(min(xx), max(xx), length.out = 2001L)
        approx(f(xs), xs, xout = c(-cut, cut), ties = mean)$y
      }, error = function(e) approx(tt, xx, xout = c(-cut, cut), ties = mean)$y)
      out[nm, ] <- inv
    }
    for (side in 1:2) if (is.na(out[nm, side]))
      msgs <- c(msgs, paste0(nm, ", ", c("lower", "upper")[side], " limit: ", .why_na(why[[side]])))
  }
  if (length(msgs)) warning("Some limits are NA:\n  ", paste(msgs, collapse = "\n  "), call. = FALSE)
  out
}

plot.profile.onls <- function(x, absVal = TRUE, level = 0.95, ...)
{
  fit <- attr(x, "original.fit")
  nm  <- names(x)
  cut <- .ci_quantile(fit, level)
  op  <- par(mfrow = grDevices::n2mfrow(length(nm))); on.exit(par(op))
  for (p in nm) {
    d  <- x[[p]]; pv <- d$par.vals[, p]
    tt <- if (absVal) abs(d$tau) else d$tau
    plot(pv, tt, type = "b", pch = 16, xlab = p, ylab = if (absVal) "|tau|" else "tau", ...)
    abline(h = if (absVal) cut else c(-cut, cut), lty = 2, col = "grey50")
    abline(v = coef(fit)[[p]], lty = 3, col = "grey50")
  }
  invisible(NULL)
}

confint.onls <- function(object, parm, level = 0.95, method = c("profile", "wald", "boot"),
                         k = 200, max_tries = 10 * k, verbose = FALSE, ...)
{
  method <- match.arg(method)
  cf <- coef(object); pn <- names(cf)
  free <- if (!is.null(object$fixed)) !as.logical(object$fixed) else rep(TRUE, length(cf))
  if (missing(parm)) parm <- pn[free]
  if (is.numeric(parm)) parm <- pn[parm]
  if (any(!parm %in% pn)) stop("Unknown parameter(s) in `parm`: ", paste(setdiff(parm, pn), collapse = ", "))
  if (any(!free[match(parm, pn)])) {
    warning("Fixed parameters have no confidence interval: ", paste(parm[!free[match(parm, pn)]], collapse = ", "))
    parm <- parm[free[match(parm, pn)]]
  }

  if (!length(parm)) return(matrix(NA_real_, 0L, 2L, dimnames = list(NULL, .ci_names(level))))
  
  if (method == "wald") {
    q  <- .ci_quantile(object, level)
    se <- object$std_errors[parm]
    out <- cbind(cf[parm] - q * se, cf[parm] + q * se)
    dimnames(out) <- list(parm, .ci_names(level))
    return(out)
  }

  if (method == "profile") {
    prof <- profile(object, which = parm, alphamax = (1 - level) / 4, ...)
    return(confint(prof, parm = parm, level = level))
  }

  ## nonparametric case bootstrap: rows resampled together with their weights and sigmas
  boot <- matrix(NA_real_, k, length(cf), dimnames = list(NULL, pn))
  need_ortho <- all(object$ortho %in% TRUE)
  n_ok <- 0L; tries <- 0L
  while (n_ok < k && tries < max_tries) {
    tries <- tries + 1L
    idx <- sample.int(object$n, replace = TRUE)
    fit <- tryCatch(.refit(object, idx = idx), error = function(e) NULL)
    if (is.null(fit) || fit$convergence != 0L) next
    if (need_ortho && !all(fit$ortho %in% TRUE)) next
    n_ok <- n_ok + 1L
    boot[n_ok, ] <- coef(fit)
    if (verbose) { cat(if (n_ok %% 10L == 0L) n_ok else "."); if (n_ok %% 50L == 0L) cat("\n"); flush.console() }
  }
  if (verbose) cat("\n")
  if (n_ok < k) warning("Only ", n_ok, " of ", k, " bootstrap fits succeeded in ", tries, " attempts.")
  if (n_ok < 2L) stop("Too few successful bootstrap fits for an interval.")
  a <- (1 - level) / 2
  out <- t(apply(boot[seq_len(n_ok), parm, drop = FALSE], 2L, quantile, probs = c(a, 1 - a), names = FALSE))
  dimnames(out) <- list(parm, .ci_names(level))
  attr(out, "n_boot") <- n_ok
  out
}
