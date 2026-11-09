\name{predict.onls}
\alias{predict.onls}
\encoding{latin1}

\title{Predictions, confidence and prediction intervals from an onls fit}

\description{
Evaluates the fitted \code{\link{onls}} model at new predictor values, optionally with standard errors and with confidence intervals for the model function or prediction intervals for a new observation in the errors-in-variables model.
}

\usage{
\method{predict}{onls}(object, newdata, se.fit = FALSE, scale = NULL, df = NULL,
        interval = c("none", "confidence", "prediction"), level = 0.95,
        sigma_x = NULL, sigma_y = NULL, ...)
}

\arguments{
  \item{object}{an object returned from \code{\link{onls}}.}
  \item{newdata}{a data frame with columns for all predictors. If missing, the observed predictors of the fit are used.}
  \item{se.fit}{logical: return the standard errors of the fitted values?}
  \item{scale}{the residual standard deviation \eqn{s}. Default: \eqn{1} with known variances, otherwise \eqn{\sqrt{\mathrm{reduced\_chisq}}}{sqrt(reduced_chisq)}. If given, the standard errors of the fitted values are rescaled to it too, as in \code{\link{predict.lm}}.}
  \item{df}{the degrees of freedom of the \eqn{t} quantiles. Default: \code{Inf} (normal quantiles) with known variances, otherwise \code{df_resid}.}
  \item{interval}{type of interval: none, a confidence interval for the model function \eqn{f(x, \theta)}, or a prediction interval for a new observation \eqn{y} at the observed (error-prone) predictor \eqn{x}.}
  \item{level}{the confidence level.}
  \item{sigma_x, sigma_y}{the standard deviations of the predictor and response errors of the new observations, in any form accepted by \code{\link{onls}}: a scalar, a length-\eqn{p} vector, a \eqn{p \times p}{p x p} covariance matrix (a square \eqn{p \times p}{p x p} matrix is always taken as a covariance, also if \code{newdata} has \eqn{p} rows; give per-row standard deviations for \eqn{p} rows as a data frame), or one row of standard deviations per row of \code{newdata}; \code{sigma_y} of length 1 or \code{nrow(newdata)}. Needed for \code{interval = "prediction"} only if the fit used observation-specific values; otherwise those of the fit are used.}
  \item{...}{further arguments, ignored.}
}

\details{
With \eqn{J} the gradient of the model function with respect to the free parameters at \code{newdata} and \eqn{V} their covariance (\code{\link{vcov}}), the standard error of the fitted value is \eqn{se_{\mathrm{fit}} = \sqrt{J V J^T}}{se_fit = sqrt(J V J')}.

A prediction interval refers to a new observation whose \emph{response and predictors} are both measured with error. Its variance is the effective variance of ODR (as in the parameter covariance, see \code{\link{onls}}) plus the variance of the fitted value:
\deqn{se^2 = se_{\mathrm{fit}}^2 + s^2\left(\sigma_{y,\mathrm{new}}^2 + g^T \Sigma_{x,\mathrm{new}}\, g\right),}{se^2 = se_fit^2 + s^2 (sigma_y,new^2 + g' Sigma_x,new g)}
where \eqn{g} is the gradient of the model function with respect to the predictors at \code{newdata}, and \eqn{\Sigma_{x,\mathrm{new}}}{Sigma_x,new} the predictor covariance of the new observation. The quantiles are taken from the normal distribution for known variances (\code{known_sigma = TRUE}) and from the \eqn{t} distribution with \code{df_resid} degrees of freedom otherwise, as in \code{\link{summary.onls}} and \code{\link{confint.onls}}.

A prediction interval at the \emph{observed} predictor is the relevant one for a future measurement pair \eqn{(x, y)}. If the predictor of the new observation is known exactly, use \code{sigma_x = 0} (\emph{i.e.} the response error only). A new observation has weight 1: if the fit used non-constant \code{weights}, a warning suggests to supply its \code{sigma_y}.
}

\value{
If \code{se.fit = FALSE} and \code{interval = "none"}, a vector of fitted values. With an interval, a matrix with columns \code{fit}, \code{lwr} and \code{upr}. If \code{se.fit = TRUE}, a list with the components \code{fit} (vector or matrix as before), \code{se.fit} (the standard errors of the fitted values), \code{df} and \code{residual.scale} (\eqn{s}).
}

\author{
Andrej-Nikolai Spiess
}

\seealso{
\code{\link{onls}}, \code{\link{confint.onls}}, \code{\link{summary.onls}}.
}

\examples{
## straight line, known error variances in x and y
set.seed(1)
xt <- 1:20; x <- xt + rnorm(20, 0, 0.5); y <- 2 + 0.5 * xt + rnorm(20, 0, 0.3)
mod <- onls(y ~ a + b * x, data.frame(x, y), start = list(a = 1, b = 1),
            sigma_x = 0.5, sigma_y = 0.3)
nd <- data.frame(x = c(5, 15))
predict(mod, nd, interval = "confidence")
predict(mod, nd, interval = "prediction")
## new observations with an error-free predictor
predict(mod, nd, interval = "prediction", sigma_x = 0)
}

\keyword{models}
\keyword{nonlinear}
