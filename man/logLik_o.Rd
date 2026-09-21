\name{logLik_o}
\alias{logLik_o}
\encoding{latin1}

\title{Log-Likelihood for the orthogonal residuals}

\description{
Returns the log-likelihood as calculated from the orthogonal residuals obtained by \code{\link{residuals_o}}, including a Jacobian/normalizing-constant correction for the response and predictor precisions (\code{weights}, \code{sigma_x}, \code{sigma_y}) used by \code{\link{onls}}, so that the resulting value (and hence \code{\link{AIC}}/\code{\link{BIC}}) is comparable across \code{\link{onls}} fits made with different weighting schemes.
}

\usage{
logLik_o(object) 
}

\arguments{
\item{object}{an object returned from \code{\link{onls}}.}
}

\value{
The log-likelihood as calculated from the orthogonal residuals, with attributes \code{"df"} (\eqn{1 + q_{free}}, the number of free model parameters plus one for the residual scale) and \code{"nobs"}/\code{"nall"} (the number of observations).
}

\details{
Let \eqn{d_i} be the weighted orthogonal distance of observation \eqn{i} at the fitted parameters (\code{\link{residuals_o}}, i.e. \code{object$resid_o}), and \eqn{N} the number of observations. A naive treatment of the \eqn{d_i} as i.i.d. \eqn{N(0,\sigma^2)} draws gives the usual concentrated Gaussian log-likelihood
\deqn{\ell_0 = -\frac{N}{2}\left(\log(2\pi) + 1 - \log(N) + \log\left(\sum_i d_i^2\right)\right),}
which is what earlier versions of \code{logLik_o} returned. However, since \code{\link{onls}} already folds the response precision \eqn{Qyy_i} and predictor precision \eqn{Qx_i} into \eqn{d_i} itself (see 'Details' in \code{\link{onls}}), \eqn{\ell_0} alone omits the normalizing-constant (Jacobian) term that these precisions contribute to the underlying Gaussian density -- exactly as plain \code{\link{logLik}}/\code{\link{AIC}} on \emph{unweighted} residuals from a \emph{weighted} \code{\link{lm}}/\code{\link{nls}} fit would, if computed without R's own \eqn{+\frac{1}{2}\sum_i \log(w_i)} correction. \code{logLik_o} therefore adds the corresponding correction,
\deqn{\ell = \ell_0 \;+\; \frac{1}{2}\sum_{i=1}^{N} \log(Qyy_i) \;+\; \frac{1}{2}\sum_{i=1}^{N} \log\left|Qx_i\right|,}
where \eqn{|Qx_i|} is the determinant of the (possibly per-observation) predictor precision matrix used by \code{\link{onls}} -- a scalar for single-predictor models, and the full \eqn{p \times p} determinant (correctly accounting for any predictor correlation) for multivariate models. For the default, fully unweighted \code{\link{onls}} fit, \eqn{Qyy_i = 1} and \eqn{Qx_i = I_p} for every \eqn{i}, so both correction terms are exactly \eqn{0} and \eqn{\ell = \ell_0}; the correction only changes the value for fits using \code{weights}, \code{sigma_x}, and/or \code{sigma_y}.

Note that \eqn{d_i} itself remains a simplification: it is the square root of a sum of \eqn{p+1} squared, precision-weighted Gaussian terms (one from the response, \eqn{p} from the predictors), not a single univariate normal draw, so treating \eqn{\sum_i d_i^2 / N} as the MLE of a common residual variance (as \eqn{\ell_0} does) is itself an approximation, inherited unchanged from the original (unweighted) formula. The correction above addresses only the missing weighting/precision Jacobian term, not this deeper simplification.
}

\note{
\code{logLik_o} has no other generic functions on top, so for calculating \code{\link{AIC}}, one has to apply 
\code{AIC(logLik_o(model))}, see 'Examples'. 
The usual \code{\link{logLik}} applies to the vertical residuals of the orthogonal model. 
Because \code{logLik_o} includes the response/predictor precision correction described in 'Details', its value (and \code{AIC}/\code{BIC} derived from it) can be validly compared between two \code{\link{onls}} fits on the same data even if they used different \code{weights}, \code{sigma_x}, or \code{sigma_y}, provided the response variable and the number of observations are the same.
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
DNase1 <- subset(DNase, Run == 1)
DNase1$density <- sapply(DNase1$density, function(x) rnorm(1, x, 0.1 * x))
mod1 <- onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
             data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1))
logLik_o(mod1)

## compare AIC of vertical versus orthogonal residuals.
AIC(mod1)
AIC(logLik_o(mod1))

## Comparing an unweighted and a weighted onls() fit on the same data:
## logLik_o() includes the precision Jacobian correction, so AIC/BIC
## from the two fits are validly comparable.
mod2 <- onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
              data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1),
              sigma_x = 0.05, sigma_y = 0.1)
AIC(logLik_o(mod1))
AIC(logLik_o(mod2))
}   

\keyword{optimize}
\keyword{models}
\keyword{nonlinear}

