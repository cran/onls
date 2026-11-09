\name{logLik_o}
\alias{logLik_o}
\encoding{latin1}

\title{Log-likelihood of the errors-in-variables model}

\description{
Returns the log-likelihood of the errors-in-variables (functional) model fitted by \code{\link{onls}}, in which both the response and the predictors are observed with Gaussian errors, with the precisions given by \code{weights}, \code{sigma_x} and \code{sigma_y}. It is computed from the weighted orthogonal distances (\code{\link{residuals_o}}) at the fitted parameters and foot points.
}

\usage{
logLik_o(object) 
}

\arguments{
\item{object}{an object returned from \code{\link{onls}}.}
}

\value{
An object of class \code{c("logLik_o", "logLik")}: the log-likelihood, with attributes \code{"df"} (the number of free model parameters, plus one if the variance scale was estimated, see 'Details') and \code{"nobs"}/\code{"nall"} (the number of observations with non-zero weight). \code{\link{AIC}} and \code{\link{BIC}} can be applied to it.
}

\details{
Each observation \eqn{i} contributes \eqn{1 + p} Gaussian error components, the response error with precision \eqn{s^{-2} Qyy_i} and the \eqn{p} predictor errors with precision matrix \eqn{s^{-2} Qx_i}, where \eqn{Qyy_i} and \eqn{Qx_i} are the precisions used by \code{\link{onls}} and \eqn{s^2} is a common variance scale. At the fitted parameters and foot points, the exponents add up to \eqn{S = \sum_i d_i^2}, the minimized objective, and with \eqn{N} observations of non-zero weight and \eqn{k = N(1 + p)} components the log-likelihood is
\deqn{\ell = -\frac{S}{2 s^2} - \frac{k}{2}\log(2\pi s^2) + \frac{1}{2}\sum_{i=1}^{N} \log Qyy_i + \frac{1}{2}\sum_{i=1}^{N} \log|Qx_i|}
With known variances (\code{known_sigma = TRUE} in \code{\link{onls}}), \eqn{s^2 = 1}, otherwise, \eqn{s^2} is replaced by its maximum-likelihood estimate \eqn{S/k}, giving \eqn{\ell = -\frac{k}{2}(\log(2\pi S/k) + 1) + \frac{1}{2}\sum\log Qyy_i + \frac{1}{2}\sum\log|Qx_i|}, and the scale counts as one estimated parameter.

\eqn{|Qx_i|} is the determinant of the predictor precision matrix (a scalar for one predictor; it accounts for predictor correlation in the multivariate case). Observations with zero weight are excluded, as in \code{\link{logLik}} for weighted fits. When \code{sigma_x} and \code{sigma_y} are scaled by a common factor and the scale is estimated, \eqn{\ell} does not change.

The foot points \eqn{\hat\xi_i} (the incidental parameters of the functional model) are maximized out but not counted in \code{"df"}: their number grows with \eqn{N}, so the usual asymptotics of the likelihood do not apply to them, and the maximum-likelihood estimate \eqn{S/k} of the scale is biased downwards (\code{reduced_chisq} = \eqn{S/\mathrm{df\_resid}} is the usual estimate). \code{logLik_o} and the \code{AIC}/\code{BIC} derived from it are therefore meant for comparing \code{onls} fits to the \emph{same} data (same response, predictors and observations) with the same \code{weights}, \code{sigma_x} and \code{sigma_y} that differ in the model function or in the fixed or free parameters; they are not comparable with \code{\link{logLik}} of an \code{\link{nls}} fit, which assumes error-free predictors.

The ratio of the predictor to the response error variance cannot be chosen by maximizing \code{logLik_o}: in the functional model, the likelihood with an estimated scale is unbounded as this ratio tends to zero (or infinity), a known property of errors-in-variables models (Solari 1969). The ratio has to be known, as in Deming regression.

Version 0.3-0 changed the definition: earlier versions treated the \eqn{d_i} as \eqn{N} (instead of \eqn{k}) Gaussian residuals with an estimated variance, also for known variances.
}

\note{
The usual \code{\link{logLik}} applies to the vertical residuals of the orthogonal model, i.e. to an ordinary least-squares likelihood at the \emph{observed} predictors.
}

\references{
The Maximum Likelihood Solution of the Problem of Estimating a Linear Functional Relationship.\cr
Solari ME.\cr
Journal of the Royal Statistical Society, Series B (1969), \bold{31}: 372-375.
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
set.seed(123)
DNase1 <- subset(DNase, Run == 1)
DNase1$density <- sapply(DNase1$density, function(x) rnorm(1, x, 0.1 * x))
mod1 <- onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
             data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1))
logLik_o(mod1)

## AIC of the errors-in-variables model (not comparable with AIC(mod1),
## the least-squares likelihood of the vertical residuals; see 'Details').
AIC(logLik_o(mod1))

## Comparing two models fitted with the same (known) error variances:
## three- and four-parameter logistic.
m3 <- onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), data = DNase1,
           start = list(Asym = 3, xmid = 0, scal = 1), sigma_x = 0.05, sigma_y = 0.1)
m4 <- onls(density ~ y0 + (Asym - y0)/(1 + exp((xmid - log(conc))/scal)), data = DNase1,
           start = list(y0 = 0, Asym = 3, xmid = 0, scal = 1), sigma_x = 0.05, sigma_y = 0.1)
c(m3 = AIC(logLik_o(m3)), m4 = AIC(logLik_o(m4)))
}   

\keyword{optimize}
\keyword{models}
\keyword{nonlinear}

