\name{x0, y0}
\alias{x0}
\alias{y0}
\encoding{latin1}

\title{x0/y0-values from orthogonal nonlinear least squares regression}

\description{
Returns the fitted foot points \eqn{\hat\xi_i} (\code{x0}) and the corresponding fitted response \eqn{f(\hat\xi_i, \hat\theta)} (\code{y0}) obtained by \code{\link{onls}} -- the points on the fitted curve/surface that the orthogonal-distance criterion projects each observation onto. See 'Details' for the criterion actually minimized, which reduces to plain Euclidean distance only in the simplest (single-predictor, unweighted) case; for multivariate and/or weighted fits, a precision-weighted quadratic form is minimized instead.
}

\usage{
x0(object) 
y0(object) 
}

\arguments{
\item{object}{an object returned from \code{\link{onls}}.}
}

\details{
For a single-predictor (\eqn{p=1}), unweighted \code{\link{onls}} fit, the foot point \eqn{\hat\xi_i} (returned by \code{x0}) is obtained by minimizing the plain Euclidean distance between the observation \eqn{(x_i, y_i)} and a point \eqn{(\xi_i, f(\xi_i,\hat\theta))} on the fitted curve,
\deqn{\min_{\xi_i} \sqrt{(x_i-\xi_i)^2 + \left[y_i-f(\xi_i,\hat\theta)\right]^2}}
More generally -- for multivariate models (\eqn{p>1}) and/or fits using \code{weights}, \code{sigma_x}, or \code{sigma_y} -- \code{\link{onls}} instead minimizes a precision-weighted quadratic form (see 'Details' in \code{\link{onls}} for the full construction),
\deqn{\min_{\xi_i} \; Qyy_i \left[y_i - f(\xi_i,\hat\theta)\right]^2 + (\xi_i-x_i)^T Qx_i (\xi_i-x_i)}
which reduces to the plain Euclidean distance above only when \eqn{p=1} and \eqn{Qyy_i = Qx_i = 1} (the default, unweighted case) -- Euclidean distance is therefore a special case of the criterion minimized, not the general rule.

\code{x0} returns \eqn{\hat\xi_i}: a length-\eqn{n} vector for single-predictor models, or an \eqn{n \times p} matrix (one column per predictor, in the order of the predictors in the formula, as in \code{object$X_obs}) for multivariate models. \code{y0} returns \eqn{f(\hat\xi_i,\hat\theta)}, always a length-\eqn{n} vector, since the response is univariate regardless of \eqn{p}.

Values are returned in the original observation order (after \code{subset}), with \code{NA} re-inserted at the positions of observations removed by \code{na.action = na.exclude}, so that they line up with the rows of the data (they are taken from the \code{xi} and \code{yi} components of the fit). Before version 0.2-1, single-predictor models returned them sorted by the predictor.
}

\value{
For single-predictor models, a length-\eqn{n} vector of \eqn{\hat\xi_i} (\code{x0}) or \eqn{f(\hat\xi_i,\hat\theta)} (\code{y0}) values. For multivariate models, \code{x0} returns an \eqn{n \times p} matrix (one column per predictor); \code{y0} still returns a length-\eqn{n} vector, since the response remains univariate.
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
DNase1 <- subset(DNase, Run == 1)
DNase1$density <- sapply(DNase1$density, function(x) rnorm(1, x, 0.1 * x))
mod <- onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
             data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1))
x0(mod)
y0(mod)
}   

\keyword{optimize}
\keyword{models}
\keyword{nonlinear}

