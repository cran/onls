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

## ---- accessors for the fitted object
.pnames <- function(object) names(object$coefficients)
.xnames <- function(object) colnames(object$X_obs)
.yname  <- function(object) as.character(object$formula[[2L]])

## ---- predictor precisions: one representation for every form of `sigma_x` ----------
## Lx is an n x p x p array of upper Cholesky factors, t(Lx[i, , ]) %*% Lx[i, , ] = Qx_i, so that the
## predictor part of the objective is |L_i delta_i|^2. The helpers below work on all observations at
## once (column-wise loops over the p predictors only).
.L_mul  <- function(L, D) {                       # rows L_i d_i
  p <- dim(L)[3L]; out <- D * 0
  for (a in seq_len(p)) out[, a] <- .rowSums(matrix(L[, a, ], ncol = p) * D, nrow(D), p)
  out
}
.Lt_mul <- function(L, V) {                       # rows L_i' v_i
  p <- dim(L)[3L]; out <- V * 0
  for (b in seq_len(p)) out[, b] <- .rowSums(matrix(L[, , b], ncol = p) * V, nrow(V), p)
  out
}
.Q_mul  <- function(L, D) .Lt_mul(L, .L_mul(L, D))           # rows Qx_i d_i
.Q_quad <- function(L, D) .rowSums(.L_mul(L, D)^2, nrow(D), ncol(D))   # d_i' Qx_i d_i
.Lt_solve <- function(L, G) {                     # rows z_i with L_i' z_i = g_i (forward substitution)
  p <- dim(L)[3L]; z <- G * 0
  for (a in seq_len(p)) {
    s <- G[, a]
    if (a > 1L) s <- s - .rowSums(matrix(L[, seq_len(a - 1L), a], ncol = a - 1L) * z[, seq_len(a - 1L), drop = FALSE],
                                  nrow(G), a - 1L)
    z[, a] <- s / L[, a, a]
  }
  z
}
.Q_diag <- function(L) {                          # n x p: diagonal of Qx_i
  p <- dim(L)[3L]; out <- matrix(0, dim(L)[1L], p)
  for (b in seq_len(p)) out[, b] <- .rowSums(matrix(L[, , b], ncol = p)^2, dim(L)[1L], p)
  out
}
.Q_logdet <- function(L) {                        # log |Qx_i|
  p <- dim(L)[3L]; 2 * .rowSums(log(matrix(sapply(seq_len(p), function(a) L[, a, a]), ncol = p)), dim(L)[1L], p)
}

## Build precisions. Returns Q_yy (length n), Lx (n x p x p), per_obs (sigma_x given per
## observation) and is_diag (no correlated predictor errors).
.build_Q <- function(weights, sigma_x, sigma_y, pred_names, n) {
  p <- length(pred_names)
  if (!is.numeric(sigma_y) || anyNA(sigma_y)) stop("`sigma_y` must be numeric without NA.")
  if (!is.null(sigma_x) && (!is.numeric(sigma_x) || anyNA(sigma_x)))
    stop("`sigma_x` must be numeric without NA.")
  
  if (length(sigma_y) == 1L) {
    if (sigma_y <= 0) stop("`sigma_y` must be a single positive number (or a length-", n, " vector).")
    Qyy_base <- rep(1 / sigma_y^2, n)
  } else if (length(sigma_y) == n) {
    if (any(sigma_y <= 0)) stop("`sigma_y` values must be positive.")
    Qyy_base <- 1 / sigma_y^2
  } else stop("`sigma_y` must be a single positive number or a length-", n, " vector.")
  Q_yy <- Qyy_base * weights      # weights act on the y-side only
  
  diagL <- function(S) {          # n x p standard deviations -> diagonal factors 1/sd
    L <- array(0, c(n, p, p))
    for (a in seq_len(p)) L[, a, a] <- 1 / S[, a]
    L
  }
  per_obs <- FALSE; is_diag <- TRUE
  if (is.null(sigma_x)) {
    L <- diagL(matrix(1, n, p))
  } else if (is.matrix(sigma_x) && all(dim(sigma_x) == c(n, p))) {           # per observation
    if (any(sigma_x <= 0)) stop("`sigma_x` values must be positive.")
    L <- diagL(sigma_x); per_obs <- TRUE
  } else if (!is.matrix(sigma_x) && length(sigma_x) == n && p == 1L && n > 1L) {   # per observation, p = 1
    if (any(sigma_x <= 0)) stop("`sigma_x` values must be positive.")
    L <- diagL(matrix(sigma_x, ncol = 1L)); per_obs <- TRUE
  } else if (is.matrix(sigma_x)) {                                             # global p x p covariance
    if (!all(dim(sigma_x) == c(p, p)))
      stop("`sigma_x` matrix must be ", p, " x ", p, " (global) or ", n, " x ", p,
           " (per-observation, diagonal). Got ", nrow(sigma_x), " x ", ncol(sigma_x), ".")
    if (!isSymmetric(unname(sigma_x))) stop("`sigma_x` matrix must be symmetric.")
    if (any(eigen(sigma_x, symmetric = TRUE, only.values = TRUE)$values <= 0))
      stop("`sigma_x` matrix must be positive definite.")
    R <- chol(solve(sigma_x))                                                  # t(R) %*% R = Qx
    L <- array(rep(R, each = n), c(n, p, p))
    is_diag <- !any(abs(R[row(R) != col(R)]) > 1e-12 * max(abs(R)))
  } else if (length(sigma_x) == p || length(sigma_x) == 1L) {                  # global diagonal / scalar
    if (any(sigma_x <= 0)) stop("`sigma_x` values must be positive.")
    L <- diagL(matrix(rep_len(sigma_x, p), n, p, byrow = TRUE))
  } else stop("`sigma_x` must be NULL, a scalar, a length-", p, " vector, a ", p, " x ", p,
              " matrix, or an ", n, " x ", p, " matrix (per-observation).")
  dimnames(L) <- list(NULL, pred_names, pred_names)
  list(Q_yy = Q_yy, Lx = L, per_obs = per_obs, is_diag = is_diag)
}

## ---- main -------------------------------------------------------------------

## extend / window (p = 1 only): OPTIONAL bounds on the foot points. Default NULL
## = unbounded, which is what a true joint ODR needs.
onls <- function(formula, data, start = NULL, weights = NULL, sigma_x = NULL,
                 sigma_y = 1, known_sigma = NULL, extend = NULL, window = NULL,
                 control = list(), lower = NULL, upper = NULL, fixed = NULL,
                 subset = NULL, na.action = NULL, trace = FALSE,
                 algorithm = c("auto", "structured", "dense")) {
  algorithm <- match.arg(algorithm)
  
  sigma_y_supplied <- !missing(sigma_y)
  if (!inherits(formula, "formula") || length(formula) != 3L) stop("`formula` must be two-sided.")
  if (!is.name(formula[[2]])) stop("The response must be a plain variable name (no transformations).")
  data_expr <- substitute(data)
  data <- as.data.frame(data)
  env_f <- environment(formula); if (is.null(env_f)) env_f <- parent.frame()
  ## as in nls()/lm(): `weights`, `subset`, `sigma_x` and `sigma_y` are evaluated in `data` first,
  ## then in the calling environment, so that columns of `data` can be named directly
  env_c     <- parent.frame()
  weights   <- eval(substitute(weights), data, env_c)
  subset    <- eval(substitute(subset),  data, env_c)
  sigma_x   <- eval(substitute(sigma_x), data, env_c)
  sigma_y   <- eval(substitute(sigma_y), data, env_c)
  
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
    if (is.logical(subset)) subset <- which(subset)                  # NA in a logical subset: dropped, as in nls()
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
  Lx       <- Q$Lx                                   # n x p x p, t(L_i) L_i = Qx_i
  is_diag  <- Q$is_diag
  
  ## ---- sort by first predictor (p = 1) -- only affects output order / optional bounds
  if (p == 1L) { ORD <- order(X0[, 1L]); unord <- order(ORD)
  } else       { ORD <- seq_len(n);      unord <- seq_len(n) }
  y_s <- y0[ORD]; X_s <- X0[ORD, , drop = FALSE]
  wts_s <- wts[ORD]; Qyy_s <- Qyy_orig[ORD]
  Lx_s  <- Lx[ORD, , , drop = FALSE]; Qd_s <- .Q_diag(Lx_s)
  
  ## ---- parameters, fixed, bounds
  theta0 <- setNames(as.numeric(unlist(start)), pnms)
  FIXSEL <- if (!is.null(fixed)) which(as.logical(fixed)) else integer(0)
  qfree_idx <- setdiff(seq_len(q), FIXSEL)
  ## all parameters fixed: only the foot points are optimized (used when profiling one-parameter models)
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
  fit0 <- if (nq == 0L) NULL else tryCatch(withCallingHandlers({
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
  
  ## ---- joint problem, common to both algorithms: unknowns (free theta, delta),
  ## residual r = (sqrt(Qyy_i) (y_i - f(x_i + delta_i)), L_i delta_i), S = r'r
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
  LOWER_th <- if (!is.null(LOWER)) LOWER[qfree_idx] else rep(-Inf, nq)
  UPPER_th <- if (!is.null(UPPER)) UPPER[qfree_idx] else rep( Inf, nq)
  has_bounds <- any(is.finite(c(LOWER_th, UPPER_th, delta_lower, delta_upper)))
  
  ## Exact model gradients (deriv(), else central differences). MINPACK's own forward-difference
  ## Jacobian uses a step proportional to |par_j|: once a foot-point correction delta_j is tiny,
  ## the step drops below the resolution of f(x + delta) and d r_y / d delta_j collapses to 0,
  ## so LM drives delta_j to exactly 0 (a non-orthogonal foot point).
  model_grad <- .make_model_grad(rhs_expr, pnms, xnms, f_eval, env_f)
  ii <- seq_len(n)
  
  ## control$outer_max = TOTAL budget of LM iterations (default 5000)
  outer_max <- if (is.null(control$outer_max)) 5000L else as.integer(control$outer_max)
  ftol <- if (is.null(control$ftol)) 1e-10 else control$ftol
  ptol <- if (is.null(control$ptol)) 1e-10 else control$ptol
  diag_th <- if (!is.null(control$diag) && length(control$diag) == q) control$diag[qfree_idx] else NULL
  
  ## "auto": structured for more than 100 observations or with parameter bounds (its active set
  ## reaches the bounded optimum, minpack.lm's bounded LM can stop before it), dense otherwise
  ## (the compiled minpack.lm code has less overhead per iteration); foot-point bounds
  ## (extend/window) are only available in the dense algorithm
  if (algorithm == "auto")
    algorithm <- if (n > 100L || any(is.finite(c(LOWER_th, UPPER_th)))) "structured" else "dense"
  if (algorithm == "structured" && any(is.finite(c(delta_lower, delta_upper)))) algorithm <- "dense"
  
  ## Both branches below deliver: theta_opt (all parameters), Xi_s (sorted foot points),
  ## conv (0 = converged), iters, stop_code, stop_msg, run_err (NULL or the error message), CTRL.
  if (algorithm == "dense") {
    ## ======================================================================================
    ## DENSE: the joint problem as one nonlinear least-squares problem for minpack.lm::nls.lm,
    ## with the full (n + np) x (q + np) Jacobian. Cost per iteration O((n + q)^3).
    ## ======================================================================================
    par0    <- c(theta_nls[qfree_idx], rep(0, n * p))
    LOWER_j <- c(LOWER_th, delta_lower)
    UPPER_j <- c(UPPER_th, delta_upper)
    idx_delta <- nq + seq_len(n * p)
    
    .joint_residual <- function(par) {
      th <- theta_nls; th[qfree_idx] <- par[seq_len(nq)]
      delta <- matrix(par[idx_delta], nrow = n, ncol = p)
      yhat  <- f_eval(th, X_s + delta)
      if (is.null(yhat) || any(!is.finite(yhat))) return(rep(1e10, n + n * p))
      ry <- sqrt(Qyy_s) * (y_s - yhat)
      rx <- as.vector(.L_mul(Lx_s, delta))
      c(ry, rx)
    }
    
    .joint_jac <- function(par) {
      th <- theta_nls; th[qfree_idx] <- par[seq_len(nq)]
      gr <- model_grad(th, X_s + matrix(par[idx_delta], nrow = n, ncol = p))
      sq <- sqrt(Qyy_s)
      J  <- matrix(0, n + n * p, nq + n * p)
      J[ii, seq_len(nq)] <- -sq * gr$Jth[, qfree_idx, drop = FALSE]
      for (a in seq_len(p)) {
        J[cbind(ii, nq + (a - 1L) * n + ii)] <- -sq * gr$Gx[, a]         # d r_y,i / d delta_ia
        for (b in seq_len(p))                                              # d r_x,ia / d delta_ib
          J[cbind(n + (a - 1L) * n + ii, nq + (b - 1L) * n + ii)] <- Lx_s[, a, b]
      }
      J
    }
    
    ## minpack.lm silently caps a single nls.lm call at 1024 iterations, so the budget is spent
    ## in restarts of <= 1024, each warm-started from the last iterate (resets the trust region).
    per_call  <- min(outer_max, 1024L)
    ctrl_args <- list(maxiter = per_call, maxfev = per_call * (length(par0) + 1L), ftol = ftol, ptol = ptol)
    ## MINPACK's adaptive column scaling is the default; fixed `diag` only if requested
    if (!is.null(diag_th)) ctrl_args$diag <- c(diag_th, rep(1, n * p))
    CTRL <- do.call(minpack.lm::nls.lm.control, ctrl_args)
    
    if (trace) cat("Optimizing orthogonal NLS (dense joint problem, minpack.lm)...\n")
    run_err <- NULL; ONLS <- NULL; iters_tot <- 0L; par_cur <- par0; ss_prev <- Inf
    repeat {
      fit <- tryCatch(suppressWarnings(
        minpack.lm::nls.lm(par = par_cur, fn = .joint_residual, jac = .joint_jac, control = CTRL,
                           lower = if (has_bounds) LOWER_j, upper = if (has_bounds) UPPER_j)),
        error = function(e) { run_err <<- conditionMessage(e); NULL })
      if (is.null(fit)) break
      ONLS <- fit; iters_tot <- iters_tot + fit$niter; par_cur <- fit$par
      if (fit$info != -1L) break                                  # stopped for a reason other than maxiter
      if (iters_tot >= outer_max) break                           # budget used up
      if (!is.finite(fit$deviance) || fit$deviance > ss_prev * (1 - 1e-10)) break   # restart made no progress
      ss_prev <- fit$deviance
      if (trace) cat("  restarting LM: ", iters_tot, " iterations, objective ", format(fit$deviance, digits = 10), "\n", sep = "")
    }
    
    if (!is.null(ONLS)) {
      run_err <- NULL                                            # a later restart may have failed
      theta_opt <- theta_nls; theta_opt[qfree_idx] <- ONLS$par[seq_len(nq)]
      Xi_s  <- X_s + matrix(ONLS$par[idx_delta], nrow = n, ncol = p)
      conv  <- if (ONLS$info %in% 1:4) 0L else 1L
      iters <- iters_tot; stop_code <- ONLS$info; stop_msg <- ONLS$message
    }
    
  } else if (algorithm == "structured") {
    ## ======================================================================================
    ## STRUCTURED: Levenberg-Marquardt on the same joint problem, using its block ("arrow")
    ## structure as ODRPACK does. Each delta_i only enters observation i, so the normal
    ## equations
    ##   [ A + mu Dth    B_1 ... B_n ] [ d_theta ]     [ g_theta ]
    ##   [ B_i'          D_i + mu Di ] [ d_delta_i ] = -[ g_delta_i ]
    ## with A = sum_i a_i a_i', a_i = -sqrt(Qyy_i) df/dtheta (row i of Jt), B_i = a_i u_i',
    ## u_i = -sqrt(Qyy_i) g_i, D_i = u_i u_i' + Qx_i, are solved by eliminating every p x p
    ## block D_i (Schur complement):
    ##   (A + mu Dth - sum_i h_i a_i a_i') d_theta = -g_theta + sum_i m_i a_i,
    ##   h_i = u_i' Dm_i^-1 u_i, m_i = u_i' Dm_i^-1 g_delta_i,
    ##   d_delta_i = -Dm_i^-1 g_delta_i - (Dm_i^-1 u_i) (a_i' d_theta),   Dm_i = D_i + mu Di.
    ## Cost per iteration O(n p^3 + n q^2 + q^3). Step control: MINPACK's trust region on the
    ## scaled step ||D dz|| (More 1978), with D the running maximum of the Jacobian column norms,
    ## or control$diag for theta (then 1 for delta).
    ## ======================================================================================
    CTRL <- list(maxiter = outer_max, ftol = ftol, ptol = ptol, diag = diag_th)
    th_bounded <- any(is.finite(c(LOWER_th, UPPER_th)))
    run_err <- NULL
    th_cur <- theta_nls; dl_cur <- matrix(0, n, p)
    ## Qx_i entries as an n x p x p array (for the block solves with correlated predictor errors)
    Qarr <- array(0, c(n, p, p))
    for (a in seq_len(p)) for (b in seq_len(p)) Qarr[, a, b] <- .rowSums(matrix(Lx_s[, , a] * Lx_s[, , b], ncol = p), n, p)
    sqQyy  <- sqrt(Qyy_s)
    
    ## solve M x = b after equilibration (diagonal scaling), so that very different parameter
    ## scales do not trip solve()'s condition check; NULL if M is singular
    esolve <- function(M, b) {
      d <- sqrt(abs(diag(M))); d[!is.finite(d) | d == 0] <- 1
      tryCatch(drop(solve(M / outer(d, d), b / d)) / d, error = function(e) NULL)
    }
    ## residuals and objective at (theta, delta)
    resid_at <- function(th, dl) {
      yhat <- suppressWarnings(f_eval(th, X_s + dl))
      if (is.null(yhat) || any(!is.finite(yhat))) return(NULL)
      ry <- sqQyy * (y_s - yhat)
      rx <- .L_mul(Lx_s, dl)                                       # row i = L_i delta_i
      list(ry = ry, rx = rx, S = sum(ry^2) + sum(rx^2))
    }
    ## M_i^-1 v_i for all rows of V, with M_i = Qx_i + diag(C_i): elementwise for diagonal Qx_i,
    ## otherwise by a Cholesky factorization of all M_i at once (vectorized over the observations)
    Msolve <- function(C, V) {
      if (is_diag) return(V / (Qd_s + C))
      Lc <- array(0, c(n, p, p))                                   # Lc[i, , ] lower triangle of M_i
      for (j in seq_len(p)) {
        sj <- Qarr[, j, j] + C[, j] - (if (j > 1L) .rowSums(Lc[, j, seq_len(j - 1L), drop = FALSE]^2, n, j - 1L) else 0)
        Lc[, j, j] <- sqrt(pmax(sj, .Machine$double.xmin))
        if (j < p) for (k in (j + 1L):p)
          Lc[, k, j] <- (Qarr[, k, j] - (if (j > 1L) rowSums(Lc[, k, seq_len(j - 1L), drop = FALSE] *
                                                           Lc[, j, seq_len(j - 1L), drop = FALSE]) else 0)) / Lc[, j, j]
      }
      z <- matrix(0, n, p); w <- matrix(0, n, p)                   # Lc Lc' w = v, row-wise
      for (j in seq_len(p))
        z[, j] <- (V[, j] - (if (j > 1L) rowSums(matrix(Lc[, j, seq_len(j - 1L)], n) * z[, seq_len(j - 1L), drop = FALSE]) else 0)) / Lc[, j, j]
      for (j in rev(seq_len(p)))
        w[, j] <- (z[, j] - (if (j < p) rowSums(matrix(Lc[, (j + 1L):p, j], n) * w[, (j + 1L):p, drop = FALSE]) else 0)) / Lc[, j, j]
      w
    }
    
    if (trace) cat("Optimizing orthogonal NLS (structured joint problem)...\n")
    cur <- resid_at(th_cur, dl_cur)
    if (is.null(cur)) {
      run_err <- "model not finite at the warm-start values"
    } else {
      sc_th <- if (nq) rep(0, nq) else numeric(0); sc_dl <- matrix(0, n, p)
      radius <- NULL; mu_prev <- NULL; it <- 0L; stop_code <- -1L
      while (it < outer_max) {
        it <- it + 1L
        gr <- model_grad(th_cur, X_s + dl_cur)
        Jt <- -sqQyy * gr$Jth[, qfree_idx, drop = FALSE]           # n x nq,  rows a_i'
        U  <- -sqQyy * gr$Gx                                       # n x p,   rows u_i'
        g_th <- if (nq) drop(crossprod(Jt, cur$ry)) else numeric(0)
        g_dl <- U * cur$ry + .Q_mul(Lx_s, dl_cur)                  # rows g_delta_i'
        ## scaling D: running maximum of the Jacobian column norms (MINPACK mode 1), or control$diag
        if (is.null(diag_th)) {
          sc_th <- pmax(sc_th, sqrt(colSums(Jt^2)))
          sc_dl <- pmax(sc_dl, sqrt(U^2 + Qd_s))
          sc_th[sc_th == 0] <- 1; sc_dl[sc_dl == 0] <- 1
        } else { sc_th <- diag_th; sc_dl[] <- 1 }
        A <- crossprod(Jt)
        znorm <- sqrt(sum((sc_th * th_cur[qfree_idx])^2) + sum((sc_dl * dl_cur)^2))
        if (is.null(radius)) radius <- if (znorm > 0) 100 * znorm else 100
        
        ## step for damping mu: Schur complement of the p x p blocks, then bounds (active set)
        ## With M_i = Qx_i + mu D_i^2 (Dm_i = M_i + u_i u_i') and s_i = 1 + u_i' M_i^-1 u_i, the
        ## Schur complement A + mu Dth - sum B_i Dm_i^-1 B_i' equals sum a_i a_i' / s_i + mu Dth:
        ## written this way, nothing nearly cancels when y is much more precise than x.
        step_mu <- function(mu) {
          C  <- mu * sc_dl^2                                         # damping part of M_i
          Wm <- list(W1 = Msolve(C, U), W2 = Msolve(C, g_dl - U * cur$ry))   # M^-1 u, M^-1 Qx delta
          s1 <- 1 + .rowSums(U * Wm$W1, n, p)
          tt <- (cur$ry - .rowSums(U * Wm$W2, n, p)) / s1
          d_th <- numeric(0); a_d <- rep(0, n)
          if (nq) {
            Sch <- crossprod(Jt, Jt / s1); diag(Sch) <- diag(Sch) + mu * sc_th^2
            rhs <- -drop(crossprod(Jt, tt))
            d_th <- esolve(Sch, rhs)
            if (is.null(d_th)) return(NULL)
            th0 <- th_cur[qfree_idx]
            if (th_bounded) {
            ## parameters at a bound whose descent direction points outward stay there
            free_j <- !((th0 <= LOWER_th & rhs <= 0) | (th0 >= UPPER_th & rhs >= 0))
            if (!all(free_j)) {
              d_th[!free_j] <- 0; Fj <- which(free_j)
              if (length(Fj)) { d_F <- esolve(Sch[Fj, Fj, drop = FALSE], rhs[Fj]); if (!is.null(d_F)) d_th[Fj] <- d_F }
            }
            ## parameters whose step leaves [lower, upper] are held at the bound and the
            ## others re-solved given that (projected Newton on the Schur system)
            for (pass in seq_len(nq)) {
              out_j <- free_j & (th0 + d_th < LOWER_th | th0 + d_th > UPPER_th)
              if (!any(out_j)) break
              d_th[out_j] <- pmin(UPPER_th, pmax(LOWER_th, th0 + d_th))[out_j] - th0[out_j]
              free_j[out_j] <- FALSE
              Fj <- which(free_j); Aj <- which(!free_j)
              if (length(Fj)) { d_F <- esolve(Sch[Fj, Fj, drop = FALSE], rhs[Fj] - drop(Sch[Fj, Aj, drop = FALSE] %*% d_th[Aj]))
                                if (!is.null(d_F)) d_th[Fj] <- d_F }
            }
            d_th <- pmin(UPPER_th, pmax(LOWER_th, th0 + d_th)) - th0
            }
            a_d  <- drop(Jt %*% d_th)
          }
          d_dl <- -Wm$W2 - Wm$W1 * (tt + a_d / s1)
          list(d_th = d_th, d_dl = d_dl, a_d = a_d, C = C, W1 = Wm$W1, s1 = s1,
               Sch = if (nq) Sch else NULL,
               dnorm = sqrt(sum((sc_th * d_th)^2) + sum((sc_dl * d_dl)^2)))
        }
        ## derivative of phi(mu) = ||D dz(mu)||: phi' = -(D^2 dz)' (J'J + mu D^2)^-1 (D^2 dz) / phi,
        ## with the same block elimination (Dm_i^-1 = M_i^-1 - W1_i W1_i' / s_i)
        dphi <- function(st) {
          v_th <- sc_th^2 * st$d_th; v_dl <- sc_dl^2 * st$d_dl
          w_th <- numeric(0); a_w <- 0
          if (nq) {
            w_th <- esolve(st$Sch, v_th - drop(crossprod(Jt, .rowSums(st$W1 * v_dl, n, p) / st$s1)))
            if (is.null(w_th)) return(NA_real_)
            a_w <- drop(Jt %*% w_th)
          }
          xx   <- v_dl - U * a_w
          Mx   <- Msolve(st$C, xx)
          w_dl <- Mx - st$W1 * (.rowSums(st$W1 * xx, n, p) / st$s1)
          -(sum(v_th * w_th) + sum(v_dl * w_dl)) / max(st$dnorm, 1e-300)
        }
        
        repeat {                                                   # trust-region loop (More 1978)
          ## LM parameter: 0 if the Gauss-Newton step fits the region, else ||D step(mu)|| ~ radius
          mu <- 0; st <- step_mu(0)
          if (is.null(st) || st$dnorm > 1.1 * radius) {
            ## as MINPACK's lmpar: bracket [lo, hi] for mu, Newton steps on 1/phi - 1/radius
            ## (Hebden), started from the Newton step at mu = 0 or from the previous mu
            gnorm <- sqrt(sum((g_th / sc_th)^2) + sum((g_dl / sc_dl)^2))
            lo <- 0; hi <- max(gnorm / radius, 1e-300)
            newton <- function(st, mu) {
              dp <- dphi(st)
              if (is.finite(dp) && dp < 0) mu - (st$dnorm - radius) * st$dnorm / (radius * dp) else NA_real_
            }
            if (!is.null(st)) { c0 <- newton(st, 0); if (is.finite(c0) && c0 > 0) lo <- min(c0, hi) }
            mu <- if (!is.null(mu_prev) && mu_prev > lo && mu_prev < hi) mu_prev else if (lo > 0) lo else hi * 1e-3
            for (k in 1:10) {
              st <- step_mu(mu)
              if (is.null(st)) { lo <- mu; mu <- sqrt(max(lo, hi * 1e-12) * hi); next }
              if (abs(st$dnorm - radius) <= 0.1 * radius) break
              if (st$dnorm > radius) lo <- mu else hi <- mu
              cand <- newton(st, mu)
              mu <- if (is.finite(cand) && cand > lo && cand < hi) cand
                    else if (lo > 0) sqrt(lo * hi) else hi / 10
            }
            if (is.null(st)) { stop_code <- 6L; break }
          }
          mu_prev <- mu
          th_new <- th_cur; th_new[qfree_idx] <- th_cur[qfree_idx] + st$d_th
          dl_new <- dl_cur + st$d_dl
          new <- resid_at(th_new, dl_new)
          ## reduction predicted by the linearized model: S - |r + J dz|^2
          lin_y <- cur$ry + st$a_d + .rowSums(U * st$d_dl, n, p)
          lin_x <- cur$rx + .L_mul(Lx_s, st$d_dl)
          pred  <- cur$S - (sum(lin_y^2) + sum(lin_x^2))
          act   <- if (is.null(new)) -Inf else cur$S - new$S
          rho   <- if (pred > 0) act / pred else -1
          ## radius update as in MINPACK's lmder
          if (rho <= 0.25) radius <- 0.5 * min(radius, 10 * st$dnorm)
          else if (mu == 0 || rho >= 0.75) radius <- 2 * st$dnorm
          if (rho >= 1e-4) {                                       # accept the step
            f_ok <- abs(act) <= ftol * cur$S && pred <= ftol * cur$S && rho <= 2
            x_ok <- radius <= ptol * max(znorm, .Machine$double.eps)
            th_cur <- th_new; dl_cur <- dl_new; cur <- new
            stop_code <- if (f_ok && x_ok) 3L else if (f_ok) 1L else if (x_ok) 2L else -1L
            break
          }
          ## rejected: the region has shrunk; stop if nothing more can be gained
          if (radius <= ptol * max(znorm, .Machine$double.eps)) { stop_code <- 2L; break }
          if (abs(pred) <= .Machine$double.eps * cur$S) { stop_code <- if (ftol > 0) 1L else 6L; break }
          if (radius < 1e-300) { stop_code <- 6L; break }
        }
        if (trace && (it %% 10L == 0L || stop_code != -1L))
          cat("  iteration ", it, ": objective ", format(cur$S, digits = 10), "\n", sep = "")
        if (stop_code != -1L) break
      }
      theta_opt <- th_cur
      Xi_s  <- X_s + dl_cur
      conv  <- if (stop_code %in% 1:4) 0L else 1L
      iters <- it
      stop_msg <- switch(as.character(stop_code),
        "1" = "Relative error in the sum of squares is at most `ftol'.",
        "2" = "Relative error between `par' and the solution is at most `ptol'.",
        "3" = "Conditions for `info = 1' and `info = 2' both hold.",
        "6" = "No further reduction of the sum of squares is possible (trust region exhausted).",
        "-1" = paste0("Number of iterations has reached `outer_max' == ", outer_max, "."))
    }
  }
  
  if (is.null(run_err)) {
    if (trace) cat(if (conv == 0L) "  Passed..." else "  Did not fully converge:", stop_msg, "\n")
    if (conv != 0L)
      warning("ODR optimisation did not converge (info = ", stop_code, ": ", stop_msg,
              "). Try other starting values or raise control$outer_max (total iteration budget, default 5000).")
  } else {
    theta_opt <- theta_nls
    Xi_s  <- X_s
    conv  <- 1L; iters <- NA_integer_; stop_code <- NA_integer_
    stop_msg <- paste("outer ODR optimisation failed to run:", run_err)
    warning("Outer ODR optimisation failed to run (", run_err, "); returning NLS estimates.")
  }
  theta_opt <- setNames(theta_opt, pnms)
  
  ## ---- final quantities (original observation order)
  Xi_fin <- Xi_s[unord, , drop = FALSE]
  
  fittedONLS <- f_eval(theta_opt, X0)
  
  DX_fin <- Xi_fin - X0
  Yi_fin <- f_eval(theta_opt, Xi_fin)                 # foot-point responses, original order
  r_fin  <- y0 - Yi_fin
  quad_x_fin <- .Q_quad(Lx, DX_fin)
  orth_d <- sqrt(pmax(Qyy_orig * r_fin^2 + quad_x_fin, 0))
  
  gr_fin <- model_grad(theta_opt, Xi_fin)
  Jth    <- gr_fin$Jth
  df_r <- max(sum(Qyy_orig > 0) - nq, 1L)          # zero-weight observations carry no information
  if (is.null(known_sigma)) known_sigma <- !is.null(sigma_x) || sigma_y_supplied
  known_sigma <- isTRUE(as.logical(known_sigma))
  if (known_sigma && is.null(sigma_x))
    warning("`sigma_x` was not supplied, but the variances are treated as known (`known_sigma = TRUE`, ",
            "e.g. because `sigma_y` was given): the predictor errors are taken to have standard deviation 1 ",
            "in the units of each predictor. Supply `sigma_x`, or set `known_sigma = FALSE`.")
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
  gSg   <- .rowSums(.Lt_solve(Lx, G)^2, n, p)                        # g_i' Qx_i^-1 g_i
  w_eff <- 1 / (1 / Qyy_orig + gSg)
  vcov_free <- if (nq == 0L) matrix(0, 0L, 0L) else tryCatch({
    Minf <- t(Jfree) %*% (Jfree * w_eff)
    dq   <- sqrt(abs(diag(Minf))); dq[!is.finite(dq) | dq == 0] <- 1
    raw  <- solve(Minf / outer(dq, dq)) / outer(dq, dq)              # equilibrated inverse
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
    stopCode = stop_code, stopMessage = stop_msg)
  
  OUT <- list(
    call = match.call(), data = data_expr, formula = formula, algorithm = algorithm,
    coefficients = theta_opt, std_errors = se, vcov = vcov_mat, parNLS = theta_nls,
    convergence = conv, convInfo = convInfo, objective = sum(orth_d^2), reduced_chisq = sig2,
    known_sigma = known_sigma, df_resid = df_r, n = n,
    X_obs = X0, y_obs = y0, xi = Xi_fin, yi = Yi_fin, fittedONLS = fittedONLS, orth_dist = orth_d,
    weights = wts, Q_yy = Qyy_orig, Lx = Lx, sigma_x = sigma_x, sigma_y = sigma_y, per_obs_x = Q$per_obs,
    fixed = fixed, lower = LOWER, upper = UPPER, control = control, extend = extend, window = window,
    na.action = na.info,
    dataClasses = vapply(data[c(resp_name, xnms)], function(v) class(v)[1L], character(1L)))
  class(OUT) <- c("onls", "nls")
  
  CHECK <- tryCatch(check_o(OUT, plot = FALSE), error = function(e) NULL)
  OUT$ortho <- if (!is.null(CHECK)) CHECK$Ortho else rep(NA, n)
  
  ## ---- foot-point bounds: the bounded LM of minpack.lm can stop at a non-stationary point
  ## (seen with `window`). Warn if a foot point is neither orthogonal nor at its bound.
  if (has_bounds && p == 1L && !is.null(extend)) {
    lo <- (X_s[, 1L] + delta_lower)[unord]; hi <- (X_s[, 1L] + delta_upper)[unord]
    tolb <- 1e-8 * max(diff(range(X0[, 1L])), 1)
    at_bound <- abs(Xi_fin[, 1L] - lo) < tolb | abs(Xi_fin[, 1L] - hi) < tolb
    bad <- which(!(OUT$ortho %in% TRUE) & !at_bound)
    if (length(bad))
      warning(length(bad), " foot point(s) are neither orthogonal nor at their `extend`/`window` bound: ",
              "the bounded solver stopped early. Try the fit without `window` (or without `extend`).")
  }
  
  OUT
}
