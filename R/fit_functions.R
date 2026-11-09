x0 <- function(object)
{
  if (!inherits(object, "onls")) stop("Input must be a result of 'onls'!")
  val <- object$xi                                  # foot points, original observation order
  if (ncol(val) == 1L) val <- as.vector(val)
  if (!is.null(object$na.action)) 
    val <- napredict(object$na.action, val)
  attr(val, "label") <- "x0 values from orthogonal model"
  val  
}

y0 <- function(object)
{
  if (!inherits(object, "onls")) stop("Input must be a result of 'onls'!")
  val <- as.vector(object$yi)                       # foot-point responses, original order
  if (!is.null(object$na.action)) 
    val <- napredict(object$na.action, val)
  attr(val, "label") <- "y0 values from orthogonal model"
  val  
}

deviance_o <- function(object)
{
  val <- sum(object$orth_dist^2)
  attr(val, "label") <- "Deviance (RSS) of orthogonal residuals from orthogonal model"
  val  
}

residuals_o <- function(object)
{
  val <- object$orth_dist
  if (!is.null(object$na.action)) 
    val <- napredict(object$na.action, val)
  attr(val, "label") <- "Orthogonal residuals from orthogonal model"
  val  
}

## Log-likelihood of the errors-in-variables model at the ODR solution. Every observation i
## contributes 1 + p Gaussian components: y_i ~ N(f(xi_i), s2 / Qyy_i), x_i ~ N(xi_i, s2 Qx_i^-1).
## With S = sum d_i^2, k = N (1 + p):
##   l = -S / (2 s2) - k/2 log(2 pi s2) + 1/2 sum log Qyy_i + 1/2 sum log|Qx_i|,
## s2 = 1 for known variances, else its ML estimate S / k (then l = -k/2 (log(2 pi S / k) + 1) + ...).
## The foot points xi_i are not counted in "df": they are the same n p nuisance parameters for
## every onls fit of the same data, so AIC/BIC compare such fits only.
logLik_o <- function(object)
{
  keep  <- object$Q_yy > 0                          # zero-weight observations are excluded
  S     <- sum(object$orth_dist[keep]^2)
  N     <- sum(keep); p <- ncol(object$X_obs); k <- N * (1 + p)
  known <- isTRUE(object$known_sigma)
  s2    <- if (known) 1 else S / k
  val   <- -S / (2 * s2) - k / 2 * log(2 * pi * s2) +
    0.5 * sum(log(object$Q_yy[keep])) + 0.5 * sum(.Q_logdet(object$Lx)[keep])
  attr(val, "df") <- .q_free(object) + (if (known) 0L else 1L)
  attr(val, "nobs") <- attr(val, "nall") <- N
  attr(val, "label") <- "Log-likelihood of the errors-in-variables model (orthogonal distances)"
  class(val) <- c("logLik_o", "logLik")
  val
}

print.logLik_o <- function(x, ...) 
{
  cat("'log Lik.' ", paste(format(c(x), digits = 7), collapse = ", "), " (df=", format(attr(x, "df")), ")\n", "\"", attr(x, "label"), "\"\n", sep = "")
  invisible(x)
}

## check_o(): is every foot point a stationary point of the ODR objective?
##
## In whitened coordinates (L_i x, sqrt(Qyy_i) y), with L_i'L_i = Qx_i, the weighted orthogonal
## distance is a Euclidean distance, and the stationarity condition
##   Qx_i (xi_i - x_i) = Qyy_i (y_i - f(xi_i)) grad_x f(xi_i)
## says that the residual vector v_i is orthogonal to the tangent space T_i of the model surface.
## alpha_i = angle between v_i and T_i: 90 degrees at a stationary foot point. For one predictor
## and unit precisions this is the classical tangent angle; the same criterion now covers weighted,
## multivariate and correlated (full sigma_x) fits.
check_o <- function(object, plot = TRUE, tol_deg = 0.05, tol_kkt = NULL) {
  if (!inherits(object, "onls")) stop("Input must be a result of 'onls'!")
  if (!is.null(tol_kkt))
    warning("`tol_kkt` is no longer used: all fits are checked by the angle `tol_deg` in whitened coordinates.")
  
  X  <- object$X_obs; Y <- object$y_obs; Xi <- object$xi
  n  <- nrow(X); p <- ncol(X); xnms <- .xnames(object)
  theta <- object$coefficients; pnms <- .pnames(object)
  env_f <- environment(object$formula); if (is.null(env_f)) env_f <- globalenv()
  fx <- function(M) .eval_nl(object$formula[[3L]], theta, pnms, M, xnms, env_f, strict = TRUE)
  Yi <- if (!is.null(object$yi)) object$yi else fx(Xi)
  G  <- .grad_x(fx, Xi)                                  # n x p slopes at the foot points
  Qyy <- if (!is.null(object$Q_yy)) object$Q_yy else rep(1, n)
  
  dlt <- Xi - X; r <- Y - Yi
  ## scale for "residual of zero length": size of the data, not their distance from the origin
  spx <- apply(X, 2L, function(v) max(abs(v), diff(range(v)), 1e-300))
  spy <- max(abs(Y), diff(range(Y)), 1e-300)
  alpha <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    if (Qyy[i] <= 0) next                                # zero weight: foot point = observation
    Q  <- crossprod(matrix(object$Lx[i, , ], p, p)); g <- G[i, ]; d <- dlt[i, ]
    k  <- Qyy[i] * r[i] * g - drop(Q %*% d)              # T'v: zero at a stationary foot point
    vv <- sum(d * drop(Q %*% d)) + Qyy[i] * r[i]^2       # |v|^2 = d_i^2
    sc <- sum(spx * drop(Q %*% spx)) + Qyy[i] * spy^2
    if (!is.finite(vv) || vv <= .Machine$double.eps * sc) next    # observation on the surface
    TT <- Q + Qyy[i] * tcrossprod(g)
    proj <- sum(k * solve(TT, k))                        # |projection of v onto T|^2
    alpha[i] <- acos(min(sqrt(max(proj, 0) / vv), 1)) * 180 / pi
  }
  ok <- is.na(alpha) | abs(alpha - 90) < tol_deg
  
  if (plot) {
    plot(alpha, type = "o", pch = 16, cex = 1.2, col = ifelse(ok, "black", "darkred"),
         xlab = "Observation", ylab = "Alpha (degrees)",
         main = "Angle between residual and model surface (90 = orthogonal)")
    abline(h = 90, lty = 2, col = "grey50")
  }
  
  nameY <- .yname(object)
  if (p == 1L) {
    OUT <- data.frame(X[, 1L], Xi[, 1L], Y, Yi, alpha, G[, 1L], ok)
    colnames(OUT) <- c(xnms, "x0", nameY, "y0", "alpha", "df/dx", "Ortho")
  } else {
    OUT <- data.frame(X, Xi, Y, Yi, alpha, G, ok)
    colnames(OUT) <- c(xnms, paste0("x0_", xnms), nameY, "y0", "alpha", paste0("df/dx_", xnms), "Ortho")
  }
  OUT
}
