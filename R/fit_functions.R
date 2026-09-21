x0 <- function(object)
{
  if (class(object)[1] != "onls") stop("Input must be a result of 'onls'!")
  val <- object$x0
  if (is.null(dim(val))) val <- as.vector(val)   # only flatten true vectors (p == 1)
  if (!is.null(object$na.action)) 
    val <- napredict(object$na.action, val)
  lab <- "x0 values from orthogonal model"
  attr(val, "label") <- lab
  val  
}

y0 <- function(object)
{
  if (class(object)[1] != "onls") stop("Input must be a result of 'onls'!")
  val <- as.vector(object$y0)   # unchanged -- y0 is always a plain vector, even for p > 1
  if (!is.null(object$na.action)) 
    val <- napredict(object$na.action, val)
  lab <- "y0 values from orthogonal model"
  attr(val, "label") <- lab
  val  
}

deviance_o <- function(object)
{
  val <- sum(object$orth_dist^2)
  lab <- "Deviance (RSS) of orthogonal residuals from orthogonal model"
  attr(val, "label") <- lab
  val  
}

residuals_o <- function(object)
{
  val <- object$orth_dist
  if (!is.null(object$na.action)) 
    val <- napredict(object$na.action, val)
  lab <- "Orthogonal residuals from orthogonal model"
  attr(val, "label") <- lab
  val  
}

logLik_o <- function(object)
{
  res_o <- object$resid_o
  N <- length(res_o)
  Qyy_i <- object$Q_yy
  Qx <- object$Q_x
  
  if (isTRUE(object$per_obs_x)) log_det_Qx_i <- rowSums(log(Qx))
  else {
    log_det_Qx <- if (is.matrix(Qx)) as.numeric(determinant(Qx, logarithm = TRUE)$modulus) else sum(log(Qx))
    log_det_Qx_i <- rep(log_det_Qx, N)
  }
  
  jacobian_corr <- 0.5 * sum(log(Qyy_i)) + 0.5 * sum(log_det_Qx_i)
  val_o <- -N * (log(2 * pi) + 1 - log(N) + log(sum(res_o^2))) / 2 + jacobian_corr
  q_free <- object$q - if (!is.null(object$fixed)) sum(as.logical(object$fixed)) else 0L
  val <- val_o
  attr(val, "df") <- 1L + q_free
  attr(val, "nobs") <- attr(val, "nall") <- N
  attr(val, "label") <- "Log-likelihood using orthogonal residuals from orthogonal model"
  class(val) <- c("logLik_o", "logLik")
  val
}

print.logLik_o <- function(x, ...) 
{
  cat("'log Lik.' ", paste(format(c(x), digits = 7), collapse = ", "), " (df=", format(attr(x, "df")), ")\n", "\"", attr(x, "label"), "\"\n", sep = "")
  invisible(x)
}

## check_o(): verify that every foot point (x0, y0) is orthogonal to the fitted model.
##
## Unweighted:  angle between the residual vector r = (x - x0, y - y0) and the model
##              tangent t = (1, m) in the (x_k, y) plane; 90 degrees = orthogonal.
##              alpha = atan2(|m*dx - dy|, |dx + m*dy|)
##              (identical to the old atan(|(m - s)/(1 + m*s)|) with s = dy/dx, but with no
##              division by dx, so dx = 0 (foot point == observation, e.g. where m = 0) no
##              longer gives 0/0 = NaN.)
## Weighted:    relative KKT residual  Qyy*(y-y0)*m  vs  Qx*(x0-x)   (0 = exact).
##
## Self-contained: slopes are computed here (no dependency on the signature of package
## helpers), evaluated in environment(formula).
check_o <- function(object, plot = TRUE, tol_deg = 0.05, tol_kkt = 1e-3) {
  if (class(object)[1] != "onls") stop("Input must be a result of 'onls'!")
  
  X     <- object$pred          # vector (p == 1) or matrix n x p (p > 1)
  Y     <- object$resp
  X0    <- x0(object)
  Y0    <- y0(object)
  theta <- object$coefficients
  pnms  <- object$param_names
  xnms  <- object$pred_names
  nameY <- attr(Y, "name")
  
  FORMULA <- object$formula[[3]]
  env_f   <- environment(object$formula); if (is.null(env_f)) env_f <- globalenv()
  
  is_weighted <- isTRUE(object$known_sigma) || !is.null(object$sigma_x) ||
    (!is.null(object$weights) && length(unique(object$weights)) > 1L)
  
  ## d f / d x_k at the foot points (central differences, step relative to |x|)
  slopes <- function(Xm) {
    fx <- function(M) {
      env <- c(setNames(as.list(theta), pnms), setNames(lapply(seq_len(ncol(M)), function(j) M[, j]), xnms))
      eval(FORMULA, envir = env, enclos = env_f)
    }
    G <- matrix(0, nrow(Xm), ncol(Xm), dimnames = list(NULL, xnms))
    for (j in seq_len(ncol(Xm))) {
      sc <- mean(abs(Xm[, j])); if (!is.finite(sc) || sc <= 0) sc <- 1
      h  <- .Machine$double.eps^(1/3) * pmax(abs(Xm[, j]), sc)
      up <- Xm; up[, j] <- Xm[, j] + h
      dn <- Xm; dn[, j] <- Xm[, j] - h
      G[, j] <- (fx(up) - fx(dn)) / (2 * h)
    }
    G
  }
  
  ## Orthogonality measure for one axis. Returns list(val, ok).
  axis_check <- function(dx, dy, m, Qyy, Qx, sc) {
    if (is_weighted) {
      lhs <- Qyy * dy * m                       # Qyy (y-y0) m
      rhs <- -Qx * dx                           # Qx (x0-x)
      val <- abs(lhs - rhs) / pmax(abs(lhs), abs(rhs), .Machine$double.eps)
      list(val = val, ok = val < tol_kkt)
    } else {
      alpha <- atan2(abs(m * dx - dy), abs(dx + m * dy)) * 180 / pi
      ## residual of (numerically) zero length: point lies on the curve -> trivially orthogonal
      zero  <- sqrt(dx^2 + dy^2) <= sqrt(.Machine$double.eps) * sc     # sc = |x| + |y| scale
      alpha[zero] <- NA_real_
      list(val = alpha, ok = ifelse(zero, TRUE, abs(alpha - 90) < tol_deg))
    }
  }
  
  ## ---------------------------------------------------------------- p == 1
  if (is.null(dim(X))) {
    nameX <- attr(X, "name")
    ## X, Y, X0, Y0 are in SORTED predictor order; Q_yy / Q_x are stored in ORIGINAL order
    ord <- order(object$X_obs[, 1L])
    m   <- slopes(matrix(X0, ncol = 1, dimnames = list(NULL, nameX)))[, 1L]
    
    Qyy_i <- if (!is.null(object$Q_yy)) object$Q_yy[ord] else rep(1, length(X))
    qx    <- object$Q_x
    Qx11  <- if (is.null(qx)) rep(1, length(X)) else
      if (is.matrix(qx) && nrow(qx) == length(X)) qx[ord, 1L] else rep(qx[[1L]], length(X))
    
    chk <- axis_check(X - X0, Y - Y0, m, Qyy_i, Qx11, abs(X) + abs(Y))
    
    if (plot) plot(as.numeric(chk$val), type = "o", col = ifelse(chk$ok, "black", "darkred"),
                   pch = 16, cex = 1.5, xlab = "x",
                   ylab = if (is_weighted) "Relative KKT residual" else "Alpha",
                   main = if (is_weighted) "Weighted fit: orthogonality via KKT residual (0 = exact)" else "")
    
    OUT <- data.frame(X, X0, Y, Y0, chk$val, m, chk$ok)
    colnames(OUT) <- c(nameX, "x0", nameY, "y0", if (is_weighted) "rel_resid" else "alpha", "df/dx", "Ortho")
    return(OUT)
  }
  
  ## ---------------------------------------------------------------- p > 1
  n <- nrow(X); p <- ncol(X)
  X0m <- if (is.matrix(X0)) X0 else matrix(X0, ncol = 1)
  colnames(X0m) <- xnms
  slopeMat <- slopes(X0m)
  Qyy_all  <- if (!is.null(object$Q_yy)) object$Q_yy else rep(1, n)
  qx       <- object$Q_x
  
  valMat <- slopeMat; valMat[] <- NA_real_
  qualMat <- matrix(NA, n, p, dimnames = list(NULL, xnms))
  for (k in xnms) {
    Qx_k <- if (is.null(qx)) rep(1, n) else
      if (is.matrix(qx) && all(dim(qx) == c(p, p))) rep(qx[k, k], n) else
        if (is.matrix(qx) && nrow(qx) == n)           qx[, k] else rep(qx[[k]], n)
    chk <- axis_check(X[, k] - X0m[, k], Y - Y0, slopeMat[, k], Qyy_all, Qx_k, abs(X[, k]) + abs(Y))
    valMat[, k] <- chk$val; qualMat[, k] <- chk$ok
  }
  QUAL <- apply(qualMat, 1, all)
  
  if (plot) {
    matplot(valMat, type = "o", pch = 16, lty = 1, cex = 1.2, col = ifelse(QUAL, "black", "darkred"),
            ylab = if (is_weighted) "Relative KKT residual" else "Alpha", xlab = "Observation index",
            main = if (is_weighted) "Weighted fit: orthogonality via KKT residual" else "Per-predictor orthogonality")
    legend("bottomright", legend = xnms, lty = seq_len(p), bty = "n")
  }
  
  OUT <- data.frame(X, X0m, Y, Y0, valMat, slopeMat, QUAL)
  colnames(OUT) <- c(xnms, paste0("x0_", xnms), nameY, "y0",
                     paste0(if (is_weighted) "rel_resid_" else "alpha_", xnms),
                     paste0("df/dx_", xnms), "Ortho")
  OUT
}
