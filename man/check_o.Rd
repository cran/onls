\name{check_o}
\alias{check_o}
\encoding{latin1}

\title{Check the fit for orthogonality of all points}

\description{
  Checks, for every observation, whether the foot point \eqn{(x_{0i}, y_{0i})} returned by \code{\link{onls}} is a genuine stationary point of the orthogonal-distance objective that \code{\link{onls}} minimizes, \emph{i.e.} whether the first-order condition \eqn{\partial S/\partial \delta_i = 0} of the joint ODRPACK-type problem (see \code{\link{onls}}) is met.
  The check is the angle between the residual vector (from the foot point to the observation) and the model curve or surface at the foot point, which is \eqn{90^{\circ}} at a stationary foot point. For weighted fits the angle is measured after each axis has been rescaled by its precision, so the same criterion covers unweighted, weighted, multivariate and correlated-error fits. See 'Details'.
}

\usage{
  check_o(object, plot = TRUE, tol_deg = 0.05, tol_kkt = NULL) 
}

\arguments{
  \item{object}{an object returned from \code{\link{onls}}.}
  \item{plot}{logical. If \code{TRUE}, the angles \eqn{\alpha_i} are plotted against the observation index, in black for orthogonal points and dark red otherwise, with a reference line at \eqn{90^{\circ}}.}
  \item{tol_deg}{tolerance in degrees: a point is orthogonal if \eqn{|\alpha_i - 90^{\circ}| < } \code{tol_deg}. Default \code{0.05}, i.e. \eqn{89.95^{\circ} < \alpha_i < 90.05^{\circ}}.}
  \item{tol_kkt}{no longer used (the relative KKT residual of version 0.2 has been replaced by the angle in whitened coordinates); supplying it gives a warning.}
}

\details{
\bold{Stationarity condition.} \code{\link{onls}} minimizes \eqn{S = \sum_i \left[ Qyy_i (y_i - f(x_i + \delta_i, \theta))^2 + \delta_i^T Qx_i \delta_i \right]} jointly over the parameters and the foot-point corrections \eqn{\delta_i = x_{0i} - x_i}. Setting \eqn{\partial S/\partial \delta_i = 0} gives \deqn{Qx_i \, (x_{0i} - x_i) = Qyy_i \, (y_i - y_{0i}) \, g_i \quad \text{with} \quad g_i = \nabla_x f(x_{0i}, \hat\theta)}

\bold{Whitened coordinates.} Let \eqn{L_i} be a factor with \eqn{L_i^T L_i = Qx_i}. In the coordinates \eqn{(L_i x, \; Qyy_i^{1/2} y)} the weighted distance \eqn{d_i} of \code{\link{onls}} is an ordinary Euclidean distance, and the condition above states that the residual vector
\deqn{v_i = \left(L_i (x_i - x_{0i}), \; Qyy_i^{1/2} (y_i - y_{0i})\right)}
is orthogonal to the tangent space of the model surface at the foot point, spanned by the columns of \eqn{T_i = (L_i^T, \; Qyy_i^{1/2} g_i)^T}. \code{check_o} computes the angle between \eqn{v_i} and this tangent space,
\deqn{\cos \alpha_i = \frac{\sqrt{k_i^T (T_i^T T_i)^{-1} k_i}}{\|v_i\|}; k_i = T_i^T v_i = Qyy_i (y_i - y_{0i}) g_i - Qx_i (x_{0i} - x_i); T_i^T T_i = Qx_i + Qyy_i g_i g_i^T,}
with \eqn{\|v_i\|^2 = d_i^2}. Hence \eqn{\alpha_i = 90^{\circ}} exactly when the stationarity condition holds. For a single predictor with unit precisions this is the classical angle between the residual vector and the tangent \eqn{(1, m_i)}, \eqn{m_i = df/dx}, as in earlier versions. The criterion is scale-free and uses the full precision matrix, so it is valid for correlated predictor errors (a global \eqn{p \times p} \code{sigma_x}) and for predictors with very small \code{sigma_x}, for which the per-axis criteria of version 0.2 reported false failures.

A point is flagged orthogonal when \eqn{|\alpha_i - 90^{\circ}| <} \code{tol_deg}. If the residual vector has (numerically) zero length, \eqn{d_i^2 \le \epsilon_{mach} (x_i^T Qx_i x_i + Qyy_i y_i^2)}, the observation lies on the fitted surface, the angle is undefined, and the point is reported with \code{alpha = NA} and \code{Ortho = TRUE}; this includes observations with zero weight.

\bold{Slopes.} The derivatives \eqn{g_i} are computed by central differences of the model formula at the foot points, evaluated in the environment of the formula, so that user-defined functions and constants used in the model are found.

\bold{Interpretation of failures.} For observations with very small residuals, the angle is sensitive to the exact position of the foot point, so that a loose convergence tolerance in \code{\link{onls}} can lead to a few flagged points although the fit is otherwise converged. Tightening \code{control = list(ftol = 1e-12, ptol = 1e-12)} in \code{\link{onls}} sharpens the angles. Foot points held at a bound set through \code{extend}/\code{window} in \code{\link{onls}} are, by construction, not stationary.
}

\value{
A data frame with one row per observation, in the original observation order (after \code{subset} and \code{na.action}):
for single-predictor models (\eqn{p=1}) the observed predictor \eqn{x_i}, the foot point \code{x0}, the observed response \eqn{y_i}, the foot-point response \code{y0}, the angle \code{alpha} in degrees (\code{NA} for an observation on the fitted curve), the model slope \code{df/dx} at the foot point, and a logical \code{Ortho};
for multivariate models (\eqn{p>1}) the observed predictors, the foot-point predictors \code{x0_<name>}, the observed response, \code{y0}, \code{alpha}, the partial derivatives \code{df/dx_<name>} at the foot point, and \code{Ortho}.
If \code{plot = TRUE}, the angles are additionally plotted.
}

\author{
  Andrej-Nikolai Spiess
}   

\examples{
## Univariate fit
set.seed(123)
x <- 1:20
y <- 10 + 3*x^2
y <- sapply(y, function(a) rnorm(1, a, 0.1 * a))
DAT <- data.frame(x, y)
mod1 <- onls(y ~ a + b * x^2, data = DAT, start = list(a = 1, b = 1))
check_o(mod1)

## Multivariate fit => one check per predictor axis
set.seed(123)
n  <- 30
x1 <- runif(n, 1, 5); x2 <- runif(n, 1, 5)
z  <- 5 + 2 * x1 + 1.5 * x2^2 + rnorm(n, 0, 2)
DAT2 <- data.frame(x1 = x1, x2 = x2, z = z)
mod2 <- onls(z ~ b1 + b2 * x1 + b3 * x2^2, data = DAT2,
             start = list(b1 = 1, b2 = 1, b3 = 1))
check_o(mod2)

## Foot point equal to the observation (model slope zero at x = 0, dx = 0):
x <- c(0, 0, 5, 7, 7.5, 10, 16, 26, 30, 34, 34.5, 100)
y <- c(1265, 1263.6, 1258, 1254, 1253, 1249.8, 1237, 1218, 1220.6, 
       1213.8, 1215.5, 1212)
DAT3 <- data.frame(x, y)
mod3 <- onls(y ~ b1 + b2 * (exp(b3 * x) - 1)^2, data = DAT3, 
             start = list(b1 = 1500, b2 = -50, b3 = -0.1))
check_o(mod3, plot = FALSE)
}   

\keyword{optimize}
\keyword{models}
\keyword{nonlinear}
