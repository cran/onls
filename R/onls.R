## =============================================================================
## onls(): orthogonal / errors-in-variables nonlinear least squares
## Joint ("combined") formulation, ODRPACK-style: theta and all foot-point
## corrections delta_i = xi_i - x_i are optimised simultaneously by LM.
##
## Objective:  sum_i  Qyy_i * (y_i - f(xi_i, theta))^2  +  delta_i' Qx_i delta_i
##
## Requires: minpack.lm.   Optional: check_o() (called if it exists).
## =============================================================================

## ---- helpers ----------------------------------------------------------------

## Evaluate the model RHS. `enclos` must be the formula's environment so that
## user-defined functions / constants are found. strict = TRUE surfaces errors
## (used once up front); otherwise errors give NULL so LM can back off.
.eval_nl <- function(expr, theta, pnms, Xi, xnms, enclos, strict = FALSE) {
  env <- c(setNames(as.list(theta), pnms),
           setNames(lapply(seq_len(ncol(Xi)), function(j) Xi[, j]), xnms))
  if (strict) return(eval(expr, envir = env, enclos = enclos))
  tryCatch(eval(expr, envir = env, enclos = enclos), error = function(e) NULL)
}

## d f / d theta   (n x q), central differences, step relative to |theta|
.grad_theta <- function(fn, theta) {
  f0 <- fn(theta)
  J  <- matrix(0, length(f0), length(theta))
  for (k in seq_along(theta)) {
    h  <- .Machine$double.eps^(1/3) * max(abs(theta[k]), 1e-3)
    up <- theta; up[k] <- theta[k] + h
    dn <- theta; dn[k] <- theta[k] - h
    J[, k] <- (fn(up) - fn(dn)) / (2 * h)
  }
  J
}

## d f / d x   (n x p), central differences, step relative to |x|
.grad_x <- function(fn, Xi) {
  G <- matrix(0, nrow(Xi), ncol(Xi), dimnames = list(NULL, colnames(Xi)))
  for (j in seq_len(ncol(Xi))) {
    sc <- mean(abs(Xi[, j])); if (!is.finite(sc) || sc <= 0) sc <- 1
    h  <- .Machine$double.eps^(1/3) * pmax(abs(Xi[, j]), sc)   # floor = data scale, not ~0
    Xup <- Xi; Xup[, j] <- Xi[, j] + h
    Xdn <- Xi; Xdn[, j] <- Xi[, j] - h
    G[, j] <- (fn(Xup) - fn(Xdn)) / (2 * h)
  }
  G
}

## Model gradients at (theta, Xi): list(Jth = n x q, Gx = n x p).
## Exact (symbolic, via deriv()) when the RHS is differentiable, else central FD.
.make_model_grad <- function(rhs_expr, pnms, xnms, f_eval, env_f) {
  fd <- function(th, Xi) list(Jth = .grad_theta(function(t) f_eval(t, Xi), th),
                              Gx  = .grad_x(function(Xm) f_eval(th, Xm), Xi))
  dexpr <- tryCatch(deriv(rhs_expr, c(pnms, xnms)), error = function(e) NULL)
  if (is.null(dexpr)) return(fd)
  function(th, Xi) {
    env <- c(setNames(as.list(th), pnms),
             setNames(lapply(seq_len(ncol(Xi)), function(j) Xi[, j]), xnms))
    g <- attr(tryCatch(eval(dexpr, envir = env, enclos = env_f), error = function(e) NULL), "gradient")
    if (is.null(g) || nrow(g) != nrow(Xi) || any(!is.finite(g))) return(fd(th, Xi))
    list(Jth = g[, pnms, drop = FALSE], Gx = g[, xnms, drop = FALSE])
  }
}

## Build precisions. Returns Q_yy (length n), Q_x (see below), Sigma_x, flags.
##   Q_x: length-p named vector (global diagonal), n x p matrix (per-obs
##        diagonal) or p x p matrix (global, correlated).
.build_Q <- function(weights, sigma_x, sigma_y, pred_names, n) {
  p <- length(pred_names)
  
  if (length(sigma_y) == 1L) {
    if (!is.numeric(sigma_y) || sigma_y <= 0)
      stop("`sigma_y` must be a single positive number (or a length-", n, " vector).")
    Qyy_base <- rep(1 / sigma_y^2, n)
  } else if (length(sigma_y) == n) {
    if (any(sigma_y <= 0)) stop("`sigma_y` values must be positive.")
    Qyy_base <- 1 / sigma_y^2
  } else stop("`sigma_y` must be a single positive number or a length-", n, " vector.")
  
  Q_yy <- Qyy_base * weights      # weights act on the y-side only
  
  if (is.null(sigma_x))
    return(list(Q_yy = Q_yy, per_obs = FALSE, is_diag = TRUE,
                Q_x = setNames(rep(1, p), pred_names),
                Sigma_x = setNames(rep(1, p), pred_names)))
  
  if (is.matrix(sigma_x) && all(dim(sigma_x) == c(n, p))) {       # per-obs, diag
    if (any(sigma_x <= 0)) stop("`sigma_x` values must be positive.")
    Sx <- sigma_x; colnames(Sx) <- pred_names
    return(list(Q_yy = Q_yy, per_obs = TRUE, is_diag = TRUE, Q_x = 1 / Sx^2, Sigma_x = Sx))
  }
  
  if (!is.matrix(sigma_x) && length(sigma_x) == n && p == 1L) {   # per-obs, p = 1
    if (any(sigma_x <= 0)) stop("`sigma_x` values must be positive.")
    Sx <- matrix(sigma_x, ncol = 1, dimnames = list(NULL, pred_names))
    return(list(Q_yy = Q_yy, per_obs = TRUE, is_diag = TRUE, Q_x = 1 / Sx^2, Sigma_x = Sx))
  }
  
  if (is.matrix(sigma_x)) {                                       # global p x p
    if (!all(dim(sigma_x) == c(p, p)))
      stop("`sigma_x` matrix must be ", p, " x ", p, " (global) or ", n, " x ", p,
           " (per-observation, diagonal). Got ", nrow(sigma_x), " x ", ncol(sigma_x), ".")
    if (!isSymmetric(unname(sigma_x))) stop("`sigma_x` matrix must be symmetric.")
    if (any(eigen(sigma_x, symmetric = TRUE, only.values = TRUE)$values <= 0))
      stop("`sigma_x` matrix must be positive definite.")
    Sigma_x <- sigma_x; dimnames(Sigma_x) <- list(pred_names, pred_names)
    offdiag <- any(abs(Sigma_x[row(Sigma_x) != col(Sigma_x)]) > 1e-12)
    ## FIX: a diagonal p x p matrix must yield a length-p *vector* Q_x,
    ## otherwise matrix(Q_x, n, p, byrow = TRUE) in onls() recycles zeros.
    Qx <- if (offdiag) solve(Sigma_x) else setNames(1 / diag(Sigma_x), pred_names)
    return(list(Q_yy = Q_yy, per_obs = FALSE, is_diag = !offdiag, Q_x = Qx, Sigma_x = Sigma_x))
  }
  
  if (length(sigma_x) == p) {                                     # global, diag
    if (any(sigma_x <= 0)) stop("`sigma_x` values must be positive.")
    return(list(Q_yy = Q_yy, per_obs = FALSE, is_diag = TRUE,
                Q_x = setNames(1 / sigma_x^2, pred_names),
                Sigma_x = setNames(sigma_x^2, pred_names)))
  }
  
  if (length(sigma_x) == 1L) {                                    # scalar
    if (sigma_x <= 0) stop("`sigma_x` must be positive.")
    return(list(Q_yy = Q_yy, per_obs = FALSE, is_diag = TRUE,
                Q_x = setNames(rep(1 / sigma_x^2, p), pred_names),
                Sigma_x = setNames(rep(sigma_x^2, p), pred_names)))
  }
  
  stop("`sigma_x` must be NULL, a scalar, a length-", p, " vector, a ", p, " x ", p,
       " matrix, or an ", n, " x ", p, " matrix (per-observation).")
}

## ---- main -------------------------------------------------------------------

## extend / window (p = 1 only): OPTIONAL bounds on the foot points. Default NULL
## = unbounded, which is what a true joint ODR needs.
onls <- function(formula, data, start = NULL, weights = NULL, sigma_x = NULL,
                 sigma_y = 1, known_sigma = NULL, extend = NULL, window = NULL,
                 control = list(), lower = NULL, upper = NULL, fixed = NULL,
                 subset = NULL, na.action = NULL, trace = FALSE) {
  
  if (!requireNamespace("minpack.lm", quietly = TRUE)) stop("Package 'minpack.lm' is required.")
  sigma_y_supplied <- !missing(sigma_y)
  if (!inherits(formula, "formula") || length(formula) != 3L) stop("`formula` must be two-sided.")
  if (!is.name(formula[[2]])) stop("The response must be a plain variable name (no transformations).")
  data_expr <- substitute(data)
  data <- as.data.frame(data)
  env_f <- environment(formula); if (is.null(env_f)) env_f <- parent.frame()
  
  resp_name <- as.character(formula[[2]])
  rhs_expr  <- formula[[3]]
  pnms      <- names(start)
  
  ## ---- call checks
  if (is.null(pnms) || length(pnms) == 0L) stop("`start` must be a named list of starting values.")
  if (!is.null(fixed) && length(fixed) != length(pnms)) stop("'fixed' must have same length as 'start'!")
  if (!resp_name %in% names(data)) stop("Response '", resp_name, "' not found in `data`.")
  xnms <- setdiff(intersect(all.vars(rhs_expr), names(data)), c(resp_name, pnms))
  if (length(xnms) == 0L) stop("No predictor variables detected.")
  p <- length(xnms); q <- length(pnms)
  
  ## row-aligned inputs (matrix n0 x p, or vector of length n0) follow row selection
  rows <- function(v, idx, n0, vec_ok) {
    if (is.null(v)) return(v)
    if (is.matrix(v) && nrow(v) == n0 && ncol(v) == p) return(v[idx, , drop = FALSE])
    if (vec_ok && !is.matrix(v) && length(v) == n0 && n0 > 1L) return(v[idx])
    v
  }
  
  ## ---- subset
  n0 <- nrow(data)
  if (!is.null(subset)) {
    if (is.logical(subset) && length(subset) != n0) stop("'subset' logical vector must have length nrow(data)")
    idx     <- seq_len(n0)[subset]
    data    <- data[idx, , drop = FALSE]
    weights <- rows(weights, idx, n0, TRUE)
    sigma_x <- rows(sigma_x, idx, n0, p == 1L)
    sigma_y <- rows(sigma_y, idx, n0, TRUE)
  }
  
  ## ---- weights
  n1  <- nrow(data)
  wts <- if (!is.null(weights)) as.numeric(weights) else rep(1, n1)
  if (length(wts) != n1) stop("'weights' must have the same length as the (subset) data")
  if (any(wts < 0 | is.na(wts))) stop("missing or negative weights not allowed")
  
  ## ---- NA handling
  na.action <- if (is.null(na.action)) {
    opt <- getOption("na.action"); if (is.null(opt)) "na.fail" else opt
  } else na.action
  check_frame <- data[, c(resp_name, xnms), drop = FALSE]
  check_frame$.weights <- wts
  na.info <- attr(match.fun(na.action)(check_frame), "na.action")
  if (!is.null(na.info)) {
    keep    <- setdiff(seq_len(n1), na.info)
    data    <- data[keep, , drop = FALSE]
    wts     <- wts[keep]
    sigma_x <- rows(sigma_x, keep, n1, p == 1L)
    sigma_y <- rows(sigma_y, keep, n1, TRUE)
  }
  
  y0 <- data[[resp_name]]; n <- length(y0)
  X0 <- as.matrix(data[xnms])
  
  f_eval <- function(th, Xi, strict = FALSE) .eval_nl(rhs_expr, th, pnms, Xi, xnms, env_f, strict)
  
  Q <- .build_Q(wts, sigma_x, sigma_y, xnms, n)
  Qyy_orig <- Q$Q_yy
  is_diag  <- Q$is_diag
  if (is_diag) {
    Qx_diag_mat <- if (Q$per_obs) Q$Q_x else
      matrix(Q$Q_x, nrow = n, ncol = p, byrow = TRUE, dimnames = list(NULL, xnms))
  } else Qx_global <- Q$Q_x
  
  ## ---- sort by first predictor (p = 1) -- only affects output order / optional bounds
  if (p == 1L) { ORD <- order(X0[, 1L]); unord <- order(ORD)
  } else       { ORD <- seq_len(n);      unord <- seq_len(n) }
  y_s <- y0[ORD]; X_s <- X0[ORD, , drop = FALSE]
  wts_s <- wts[ORD]; Qyy_s <- Qyy_orig[ORD]
  if (is_diag) Qx_diag_mat_s <- Qx_diag_mat[ORD, , drop = FALSE]
  
  ## ---- parameters, fixed, bounds
  theta0 <- setNames(as.numeric(unlist(start)), pnms)
  FIXSEL <- if (!is.null(fixed)) which(as.logical(fixed)) else integer(0)
  qfree_idx <- setdiff(seq_len(q), FIXSEL)
  if (length(qfree_idx) == 0L) stop("All parameters are fixed; nothing to estimate.")
  nq <- length(qfree_idx)
  
  LOWER <- if (!is.null(lower)) lower else control$lower
  UPPER <- if (!is.null(upper)) upper else control$upper
  if (!is.null(LOWER) && length(LOWER) != q) stop("`lower` must have length ", q)
  if (!is.null(UPPER) && length(UPPER) != q) stop("`upper` must have length ", q)
  
  ## fail loudly (once) if the model cannot be evaluated at the start values
  chk <- f_eval(theta0, X_s, strict = TRUE)
  if (!is.numeric(chk) || any(!is.finite(chk))) stop("Model gives non-finite values at `start`.")
  
  ## ---- NLS warm start (fixed parameters substituted into the formula)
  theta_nls <- theta0
  nls_df <- as.data.frame(X_s); nls_df[[resp_name]] <- y_s
  f_ws <- formula
  if (length(FIXSEL)) {
    rhs_ws <- do.call(substitute, list(rhs_expr, as.list(theta0[FIXSEL])))
    f_ws   <- as.formula(call("~", formula[[2]], rhs_ws), env = env_f)
  }
  ctrl_ws <- minpack.lm::nls.lm.control(maxiter = 1000L)
  if (trace) cat("Obtaining starting parameters from ordinary NLS...\n")
  fit0 <- tryCatch(withCallingHandlers({
    ## do.call() embeds the values in the call, so nlsLM needs no lookup of `nls_df`/`wts_s`
    ## (weights = nls_df$col fails: model.frame() evaluates it in environment(formula))
    args_ws <- list(formula = f_ws, data = nls_df, start = as.list(theta0[qfree_idx]), control = ctrl_ws)
    if (!is.null(weights)) args_ws$weights <- wts_s
    do.call(minpack.lm::nlsLM, args_ws)
  }, warning = function(w) {
    if (grepl("minFactor|singular gradient|NaNs produced", conditionMessage(w)))
      stop("NLS warm-start did not converge cleanly: ", conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    warning("NLS warm-start failed (", conditionMessage(e), "); using raw start values instead.")
    NULL
  })
  if (!is.null(fit0)) { cf <- coef(fit0); theta_nls[names(cf)] <- cf; if (trace) cat("  Passed.\n") }
  if (!is.null(LOWER)) theta_nls <- pmax(theta_nls, LOWER)
  if (!is.null(UPPER)) theta_nls <- pmin(theta_nls, UPPER)
  
  ## ---- joint problem: par = (free theta, vec(delta))
  Lx <- if (!is_diag) chol(Qx_global)      # t(Lx) %*% Lx == Qx_global
  
  delta_lower <- rep(-Inf, n * p); delta_upper <- rep(Inf, n * p)
  if (p == 1L && !is.null(extend)) {       # optional foot-point bounds
    rng  <- diff(range(X_s[, 1L]))
    XLOW <- min(X_s[, 1L]) - extend[1L] * rng
    XUPP <- max(X_s[, 1L]) + extend[length(extend)] * rng
    i <- seq_len(n)
    if (is.null(window) || n <= 25L) { lo <- rep(XLOW, n); hi <- rep(XUPP, n)
    } else {
      lo <- ifelse(i <= window,         XLOW, X_s[pmax(i - window + 1L, 1L), 1L])
      hi <- ifelse(i >= n - window + 1L, XUPP, X_s[pmin(i + window - 1L, n),  1L])
    }
    delta_lower <- lo - X_s[, 1L]; delta_upper <- hi - X_s[, 1L]
  }
  
  par0    <- c(theta_nls[qfree_idx], rep(0, n * p))
  LOWER_j <- c(if (!is.null(LOWER)) LOWER[qfree_idx] else rep(-Inf, nq), delta_lower)
  UPPER_j <- c(if (!is.null(UPPER)) UPPER[qfree_idx] else rep( Inf, nq), delta_upper)
  has_bounds <- any(is.finite(c(LOWER_j, UPPER_j)))
  idx_delta  <- nq + seq_len(n * p)
  
  .joint_residual <- function(par) {
    th <- theta_nls; th[qfree_idx] <- par[seq_len(nq)]
    delta <- matrix(par[idx_delta], nrow = n, ncol = p)
    yhat  <- f_eval(th, X_s + delta)
    if (is.null(yhat) || any(!is.finite(yhat))) return(rep(1e10, n + n * p))
    ry <- sqrt(Qyy_s) * (y_s - yhat)
    rx <- if (is_diag) as.vector(delta * sqrt(Qx_diag_mat_s)) else as.vector(delta %*% t(Lx))
    c(ry, rx)
  }
  
  ## Jacobian of the joint residual. MINPACK's own forward-difference Jacobian uses a
  ## step proportional to |par_j|: once a foot-point correction delta_j is tiny, the step
  ## drops below the resolution of f(x + delta) and d r_y / d delta_j collapses to 0, so LM
  ## drives delta_j to exactly 0 (a non-orthogonal foot point). Supply it explicitly.
  model_grad <- .make_model_grad(rhs_expr, pnms, xnms, f_eval, env_f)
  ii <- seq_len(n)
  .joint_jac <- function(par) {
    th <- theta_nls; th[qfree_idx] <- par[seq_len(nq)]
    gr <- model_grad(th, X_s + matrix(par[idx_delta], nrow = n, ncol = p))
    sq <- sqrt(Qyy_s)
    J  <- matrix(0, n + n * p, nq + n * p)
    J[ii, seq_len(nq)] <- -sq * gr$Jth[, qfree_idx, drop = FALSE]
    for (a in seq_len(p)) {
      J[cbind(ii, nq + (a - 1L) * n + ii)] <- -sq * gr$Gx[, a]         # d r_y,i / d delta_ia
      if (is_diag) {
        J[cbind(n + (a - 1L) * n + ii, nq + (a - 1L) * n + ii)] <- sqrt(Qx_diag_mat_s[, a])
      } else for (b in seq_len(p))
        J[cbind(n + (a - 1L) * n + ii, nq + (b - 1L) * n + ii)] <- Lx[a, b]
    }
    J
  }
  
  ## control$outer_max = TOTAL LM iteration budget (default 5000). minpack.lm silently caps a
  ## single nls.lm call at 1024 iterations, so the budget is spent in restarts of <= 1024,
  ## each warm-started from the last iterate (this also resets LM's trust region).
  outer_max <- if (is.null(control$outer_max)) 5000L else as.integer(control$outer_max)
  per_call  <- min(outer_max, 1024L)
  ctrl_args <- list(maxiter = per_call, maxfev = per_call * (length(par0) + 1L),
                    ftol = if (is.null(control$ftol)) 1e-10 else control$ftol,
                    ptol = if (is.null(control$ptol)) 1e-10 else control$ptol)
  ## MINPACK's adaptive column scaling is the default; fixed `diag` only if requested
  if (!is.null(control$diag) && length(control$diag) == q)
    ctrl_args$diag <- c(control$diag[qfree_idx], rep(1, n * p))
  CTRL <- do.call(minpack.lm::nls.lm.control, ctrl_args)
  
  if (trace) cat("Optimizing orthogonal NLS (joint theta + foot points)...\n")
  ONLS_err <- NULL; ONLS <- NULL; iters_tot <- 0L; par_cur <- par0; ss_prev <- Inf
  repeat {
    fit <- tryCatch(suppressWarnings(
      minpack.lm::nls.lm(par = par_cur, fn = .joint_residual, jac = .joint_jac, control = CTRL,
                         lower = if (has_bounds) LOWER_j, upper = if (has_bounds) UPPER_j)),
      error = function(e) { ONLS_err <<- conditionMessage(e); NULL })
    if (is.null(fit)) break
    ONLS <- fit; iters_tot <- iters_tot + fit$niter; par_cur <- fit$par
    if (fit$info != -1L) break                                  # stopped for a reason other than maxiter
    if (iters_tot >= outer_max) break                           # budget used up
    if (!is.finite(fit$deviance) || fit$deviance > ss_prev * (1 - 1e-10)) break   # restart made no progress
    ss_prev <- fit$deviance
    if (trace) cat("  restarting LM: ", iters_tot, " iterations, objective ", format(fit$deviance, digits = 10), "\n", sep = "")
  }
  
  ## ---- FIX: braces were mis-nested in the original (see review notes)
  if (!is.null(ONLS)) {
    theta_opt <- theta_nls; theta_opt[qfree_idx] <- ONLS$par[seq_len(nq)]
    Xi_s  <- X_s + matrix(ONLS$par[idx_delta], nrow = n, ncol = p)
    conv  <- if (ONLS$info %in% 1:4) 0L else 1L
    iters <- iters_tot
    if (trace) cat(if (conv == 0L) "  Passed..." else "  Did not fully converge:", ONLS$message, "\n")
  } else {
    theta_opt <- theta_nls
    Xi_s  <- X_s
    conv  <- 1L; iters <- NA_integer_
    warning("Outer ODR optimisation failed to run (", ONLS_err, "); returning NLS estimates.")
  }
  if (conv != 0L && !is.null(ONLS))
    warning("ODR optimisation did not converge (info = ", ONLS$info, ": ", ONLS$message,
            "). Try other starting values or raise control$outer_max (total iteration budget, default 5000).")
  theta_opt <- setNames(theta_opt, pnms)
  
  ## ---- final quantities (original observation order)
  Xi_fin <- Xi_s[unord, , drop = FALSE]
  
  fittedONLS <- f_eval(theta_opt, X0)
  fittedNLS  <- f_eval(theta_nls, X0)
  residONLS  <- y0 - fittedONLS
  residNLS   <- y0 - fittedNLS
  
  x0_sorted <- if (p == 1L) Xi_s[, 1L] else Xi_s
  y0_sorted <- f_eval(theta_opt, Xi_s)
  r_s   <- y_s - y0_sorted
  dx_s  <- Xi_s - X_s
  quad_x_s <- if (is_diag) rowSums(dx_s^2 * Qx_diag_mat_s) else rowSums(dx_s * (dx_s %*% Qx_global))
  resid_o_val <- sqrt(pmax(Qyy_s * r_s^2 + quad_x_s, 0))
  
  DX_fin <- Xi_fin - X0
  r_fin  <- y0 - f_eval(theta_opt, Xi_fin)
  quad_x_fin <- if (is_diag) rowSums(DX_fin^2 * Qx_diag_mat) else rowSums(DX_fin * (DX_fin %*% Qx_global))
  orth_d <- sqrt(pmax(Qyy_orig * r_fin^2 + quad_x_fin, 0))
  
  gr_fin <- model_grad(theta_opt, Xi_fin)
  Jth    <- gr_fin$Jth
  df_r <- max(n - nq, 1L)
  if (is.null(known_sigma)) known_sigma <- !is.null(sigma_x) || sigma_y_supplied
  sig2 <- sum(Qyy_orig * r_fin^2 + quad_x_fin) / df_r
  
  ## ---- degeneracy check: parameters the fitted model does not respond to at all
  ## (e.g. exp(b2 - b3*x) underflowed, so the model collapsed to a constant). Dimensionless
  ## sensitivity ||df/dtheta_j|| * max(|theta_j|, 1) / ||fitted||; LM can "converge" there.
  sens <- sqrt(colSums(gr_fin$Jth[, qfree_idx, drop = FALSE]^2)) *
    pmax(abs(theta_opt[qfree_idx]), 1) / max(sqrt(sum(fittedONLS^2)), .Machine$double.eps)
  dead <- pnms[qfree_idx][!is.finite(sens) | sens < 1e-8]
  if (length(dead))
    warning("Parameter(s) ", paste(dead, collapse = ", "), " have no measurable influence on the fitted ",
            "model: the solution is a degenerate stationary point (e.g. a model collapsed to a constant). ",
            "Results are not meaningful -- try different `start` values.")
  
  ## ---- covariance: effective-variance approximation, now valid for any p
  vcov_mat <- matrix(0, q, q, dimnames = list(pnms, pnms))
  Jfree <- Jth[, qfree_idx, drop = FALSE]
  G     <- gr_fin$Gx                                                    # n x p
  gSg   <- if (is_diag) rowSums(G^2 / Qx_diag_mat) else rowSums((G %*% Q$Sigma_x) * G)
  w_eff <- 1 / (1 / Qyy_orig + gSg)
  vcov_free <- tryCatch({
    raw <- solve(t(Jfree) %*% (Jfree * w_eff))
    if (known_sigma) raw else sig2 * raw
  }, error = function(e) {
    warning("Covariance matrix singular; SEs unavailable.")
    matrix(NA_real_, nq, nq)
  })
  vcov_mat[qfree_idx, qfree_idx] <- vcov_free
  se <- setNames(sqrt(pmax(diag(vcov_mat), 0)), pnms)
  if (length(FIXSEL)) se[FIXSEL] <- 0
  
  convInfo <- list(
    isConv = (conv == 0L), finIter = iters, finTol = CTRL$ftol,
    stopCode = if (!is.null(ONLS)) ONLS$info else NA_integer_,
    stopMessage = if (!is.null(ONLS)) ONLS$message else "outer ODR optimisation failed to run")
  
  mf <- as.data.frame(data)[, c(resp_name, xnms), drop = FALSE]
  for (nm in pnms) mf[[nm]] <- theta_opt[[nm]]
  
  OUT <- list(
    data = data_expr, call = match.call(), convInfo = convInfo, na.action = na.info,
    dataClasses = vapply(mf, function(v) class(v)[1L], character(1L)),
    model = mf, formula = formula,
    parNLS = theta_nls, parONLS = theta_opt,
    x0 = x0_sorted, y0 = y0_sorted,
    fittedONLS = fittedONLS, fittedNLS = fittedNLS,
    residONLS = residONLS, residNLS = residNLS, resid_o = resid_o_val,
    pred = if (p == 1L) X_s[, 1L] else X_s, resp = y_s,
    grad = Jth, QR = qr(Jth), weights = wts, control = CTRL,
    coefficients = theta_opt, std_errors = se, vcov = vcov_mat,
    xi = Xi_fin, x_corrections = DX_fin, orth_dist = orth_d, objective = sum(orth_d^2),
    Q_x = Q$Q_x, Sigma_x = Q$Sigma_x, Q_yy = Qyy_orig,
    is_diag_x = is_diag, per_obs_x = Q$per_obs,
    sigma_x = sigma_x, sigma_y = sigma_y, known_sigma = known_sigma,
    reduced_chisq = sig2, pred_names = xnms, param_names = pnms,
    start = start, fixed = fixed, lower = LOWER, upper = UPPER,
    y_obs = y0, X_obs = X0, resp_name = resp_name,
    n = n, p = p, q = q, df_resid = df_r, convergence = conv, iters = iters,
    type = "nonlinear")
  
  attr(OUT$pred, "name") <- if (p == 1L) xnms else paste(xnms, collapse = ", ")
  attr(OUT$resp, "name") <- resp_name
  class(OUT) <- c("onls", "nls")
  
  CHECK <- if (exists("check_o", mode = "function"))
    tryCatch(check_o(OUT, plot = FALSE), error = function(e) NULL) else NULL
  OUT$ortho <- if (!is.null(CHECK)) CHECK$Ortho else rep(NA, n)
  
  OUT
}
