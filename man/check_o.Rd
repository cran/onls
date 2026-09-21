\name{check_o}
\alias{check_o}
\encoding{latin1}

\title{Check the fit for orthogonality of all points}

\description{
  Checks, for every observation, whether the foot point \eqn{(x_{0i}, y_{0i})} returned by \code{\link{onls}} is a genuine stationary point of the orthogonal-distance objective that \code{\link{onls}} minimizes, \emph{i.e.} whether the first-order condition \eqn{\partial S/\partial \delta_i = 0} of the joint ODRPACK-type problem (see \code{\link{onls}}) is met.
  For unweighted (isotropic) fits, orthogonality is checked via the angle between the tangent to the model at \eqn{(x_{0i}, y_{0i})} and the Euclidean vector from the foot point to \eqn{(x_i, y_i)}, which should be \eqn{90^{\circ}}. 
  For weighted/heteroscedastic fits (see \code{sigma_x}, \code{sigma_y}, \code{weights} in \code{\link{onls}}), a plain right angle is no longer the correct criterion, so \code{check_o} instead checks the underlying first-order KKT stationarity condition directly. 
  Multivariate models (more than one predictor) are checked one predictor axis at a time. See 'Details'.
}

\usage{
  check_o(object, plot = TRUE, tol_deg = 0.05, tol_kkt = 0.001) 
}

\arguments{
  \item{object}{an object returned from \code{\link{onls}}.}
  \item{plot}{logical. If \code{TRUE}, the orthogonality diagnostic (\eqn{\alpha}-values for unweighted fits, or relative KKT residuals for weighted fits) is plotted for a quick overview of all points. For multivariate models, one panel per predictor is drawn.}
  \item{tol_deg}{tolerance in degrees for the unweighted criterion: a point is orthogonal if \eqn{|\alpha_i - 90^{\circ}| < } \code{tol_deg}. Default \code{0.05}, i.e. \eqn{89.95^{\circ} < \alpha_i < 90.05^{\circ}}.}
  \item{tol_kkt}{tolerance for the weighted criterion: a point is orthogonal if its relative KKT residual is smaller than \code{tol_kkt}. Default \code{0.001}.}
}

\details{
\bold{Which criterion is used.} \code{check_o} automatically chooses between two criteria, based on whether the fitted \code{\link{onls}} model is effectively weighted, i.e. if \emph{any} of the following holds:
\itemize{
  \item \code{object$known_sigma} is \code{TRUE},
  \item \code{object$sigma_x} was supplied to \code{\link{onls}},
  \item \code{object$weights} has more than one distinct value.
}
Otherwise (default \code{onls} call, unit \code{sigma_y}, no \code{sigma_x}, constant or absent \code{weights}) the fit is treated as unweighted/isotropic. Note that explicitly passing \code{sigma_y} to \code{\link{onls}} (even \code{sigma_y = 1}) sets \code{known_sigma = TRUE} and therefore selects the weighted (KKT) criterion.

\bold{Background.} \code{\link{onls}} minimizes \eqn{S = \sum_i \left[ Qyy_i (y_i - f(x_i + \delta_i, \theta))^2 + \delta_i^T Qx_i \delta_i \right]} jointly over the parameters and the foot-point corrections \eqn{\delta_i = x_{0i} - x_i}. Setting \eqn{\partial S/\partial \delta_i = 0} gives the stationarity condition \deqn{Qx_i \, (x_{0i} - x_i) = Qyy_i \, (y_i - y_{0i}) \, \nabla_x f(x_{0i}, \hat\theta),} which both criteria below test. For unit precisions it is exactly the statement that the vector from the foot point to the observation is orthogonal to the model surface.
  
\bold{Unweighted case: tangent-angle criterion.}\cr
Let \eqn{dx_i = x_i - x_{0i}}, \eqn{dy_i = y_i - y_{0i}} be the residual vector from the foot point to the observation, and \eqn{m_i = df(x, \hat\theta)/dx} evaluated at \eqn{x = x_{0i}} the slope of the tangent, whose direction vector is \eqn{(1, m_i)}. The function calculates the angle \eqn{\alpha_i} between the residual vector and the tangent,
\deqn{\alpha_i[^{\circ}] = \mathrm{atan2}\left(\left|m_i \, dx_i - dy_i\right|, \; \left|dx_i + m_i \, dy_i\right|\right) \cdot \frac{180}{\pi},}
which lies in \eqn{[0^{\circ}, 90^{\circ}]} and equals \eqn{90^{\circ}} exactly when \eqn{dx_i + m_i \, dy_i = 0}, the stationarity condition above with unit precisions. This is algebraically identical to the classical expression \eqn{\tan(\alpha_i) = |(m_i - n_i)/(1 + m_i n_i)|} with the slope \eqn{n_i = dy_i/dx_i} of the residual vector, but avoids the division by \eqn{dx_i}, so that points whose foot point coincides with the observation (\eqn{dx_i = 0}, for example where the model slope is zero, or observations at the boundary of a flat model) are evaluated correctly instead of returning \code{NaN}.
A point is flagged orthogonal when \eqn{|\alpha_i - 90^{\circ}| <} \code{tol_deg}. If the residual vector has (numerically) zero length, \eqn{\sqrt{dx_i^2 + dy_i^2} \le \sqrt{\epsilon_{mach}}\,(|x_i| + |y_i|)}, the observation lies on the fitted curve, the angle is undefined, and the point is reported with \code{alpha = NA} and \code{Ortho = TRUE}. This criterion is the appropriate one only when the underlying objective is the plain (unweighted) Euclidean distance, i.e. \eqn{Qyy_i = Qx_i = 1} in the notation of \code{\link{onls}}.
  
\bold{Weighted case: KKT-residual criterion.}\cr
Once \eqn{Qyy_i} and/or \eqn{Qx_i} rescale the response/predictor axes anisotropically (see \code{\link{onls}}), a plain right angle is no longer the correct orthogonality condition. Instead \code{check_o} checks directly the first-order stationarity condition that the foot point \eqn{x_{0i}} satisfies at the \code{\link{onls}} optimum,
\deqn{Qyy_i\,(y_i - y_{0i})\,m_i = Qx_i\,(x_{0i} - x_i),} by computing the relative residual of this equality,
\deqn{\text{rel\_resid}_i = \frac{\left|\,Qyy_i(y_i-y_{0i})m_i \,-\, Qx_i(x_{0i}-x_i)\,\right|}{\max\left(\left|Qyy_i(y_i-y_{0i})m_i\right|,\ \left|Qx_i(x_{0i}-x_i)\right|,\ \epsilon_{mach}\right)},}
where \eqn{\epsilon_{mach}} is \code{.Machine$double.eps} (used only to avoid division by zero when both terms vanish).
A point is flagged orthogonal when \eqn{\text{rel\_resid}_i <} \code{tol_kkt}. \eqn{Qyy_i} is taken from \code{object$Q_yy}; \eqn{Qx_i} is taken from \code{object$Q_x} (its diagonal, or observation-specific row, as applicable, see \code{\link{onls}}). For single-predictor models, \eqn{Qyy_i}/\eqn{Qx_i} (stored in the original observation order) are internally realigned to match the sorted-predictor order used elsewhere in the returned data frame (see 'Value'), so they are correctly paired with the corresponding observation.
  
\bold{Multivariate models (\eqn{p>1}).}\cr
One check (angle-based or KKT-residual-based, per the rule above) is performed \emph{per predictor axis} \eqn{k}, using the partial derivative \eqn{\partial f/\partial x_k} at the foot point in place of \eqn{m_i}, and (for the weighted case) the \eqn{k}-th diagonal precision element in place of \eqn{Qx_i}. For a global, correlated (non-diagonal) \code{sigma_x} covariance matrix, only the diagonal entries of the precision matrix are used for this per-axis check; the true joint stationarity condition couples all predictor axes simultaneously, so this diagnostic is an approximation in that case. The overall \code{Ortho} flag for an observation is \code{TRUE} only if every predictor axis passes its individual check.

\bold{Slopes.} The slopes \eqn{m_i} (or partial derivatives) are computed inside \code{check_o} by central differences of the model formula at the foot points, evaluated in the environment of the formula, so that user-defined functions and constants used in the model are found.

\bold{Interpretation of failures.} For observations with very small residuals, the angle \eqn{\alpha_i} is highly sensitive to the exact position of the foot point (an error in the foot point of relative size \eqn{\varepsilon} produces an angle error of roughly \eqn{\varepsilon} divided by the length of the residual vector), so that a loose convergence tolerance in \code{\link{onls}} can lead to a few flagged points although the fit is otherwise converged. Tightening \code{control = list(ftol = 1e-12, ptol = 1e-12)} in \code{\link{onls}} sharpens the angles. Foot points held at a bound set through \code{extend}/\code{window} in \code{\link{onls}} are, by construction, not stationary.
}

\value{
A data frame, with columns depending on whether the model has a single predictor or several, and on whether the fit is unweighted or weighted (see 'Details' for when each regime applies):
For single-predictor models (\eqn{p=1}): the observed predictor \eqn{x_i}, the foot point \eqn{x_{0i}}, the observed response \eqn{y_i}, the fitted foot-point response \eqn{y_{0i}}, either \code{alpha} (unweighted; \eqn{\alpha_i} in degrees, \code{NA} for an observation that lies on the fitted curve) or \code{rel_resid} (weighted; the relative KKT residual), the model slope \code{df/dx} at the foot point, and a logical \code{Ortho} that is \code{TRUE} when \eqn{|\alpha_i - 90^{\circ}| <} \code{tol_deg} (unweighted; also \code{TRUE} if \code{alpha} is \code{NA}) or when the relative KKT residual is \eqn{<} \code{tol_kkt} (weighted).
For multivariate models (\eqn{p>1}): the observed predictors \eqn{x_{i1},\dots,x_{ip}}, the foot-point predictors \code{x0_<name>} for each predictor, the observed response, the fitted foot-point response \code{y0}, either \code{alpha_<name>} or \code{rel_resid_<name>} for each predictor axis, the model partial derivative \code{df/dx_<name>} at the foot point for each predictor axis, and an overall logical \code{Ortho} that is \code{TRUE} only if the per-axis criterion passes on \emph{every} predictor axis.
If \code{plot = TRUE}, the diagnostic is additionally plotted, in black where \code{Ortho = TRUE} and dark red otherwise.

For single-predictor models, rows are returned in the internally-used sorted-predictor order (matching \code{object$pred}/\code{object$resp}/\code{x0}/\code{y0}), NOT the original row order of the input data; for multivariate models, sorting is a no-op and the original observation order is used.
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
