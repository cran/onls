\name{anova.onls}
\alias{anova.onls}
\encoding{UTF-8}

\title{Comparison of nested orthogonal regression models}

\description{
Compares nested models fitted by \code{\link{onls}} to the same data with the same measurement errors, by the change in the minimized orthogonal objective \eqn{S}: an F test when the error variances are estimated, a likelihood-ratio chi-square test when they are known.
}

\usage{
\method{anova}{onls}(object, ...)
}

\arguments{
  \item{object}{an object returned from \code{\link{onls}}.}
  \item{...}{further \code{onls} fits to be compared with \code{object}, in order of increasing (or decreasing) complexity.}
}

\details{
\code{onls} estimates the parameters by minimizing the weighted sum of squared orthogonal distances
\deqn{S = \sum_{i=1}^{n} d_i^2 = \sum_{i=1}^{n} \left[ Qyy_i \left(y_i - f(\hat\xi_i, \hat\theta)\right)^2 + (\hat\xi_i - x_i)^T Qx_i (\hat\xi_i - x_i) \right]}
(\code{objective} of the fit, see \code{\link{onls}}). The comparison of nested models is based on this criterion. As in \code{anova.nls}, each model is compared with the preceding one in the list. For a model \eqn{a} with \eqn{df_a} residual degrees of freedom and a model \eqn{b} with more parameters (\eqn{df_b < df_a}):

\emph{Estimated variances} (\code{known_sigma = FALSE}): the extra-sum-of-squares F test
\deqn{F = \frac{(S_a - S_b)/(df_a - df_b)}{S_b/df_b} \sim F_{df_a - df_b,\; df_b}}

\emph{Known variances} (\code{known_sigma = TRUE}): with the errors on their absolute scale, \eqn{S} is \eqn{-2} times the log-likelihood of the errors-in-variables model up to a constant (see \code{\link{logLik_o}}), so that \eqn{S_a - S_b} is the likelihood-ratio statistic, referred to a \eqn{\chi^2} distribution with \eqn{df_a - df_b} degrees of freedom.
The models must have been fitted to the same observations, with the same \code{weights}, \code{sigma_y} and \code{sigma_x} for the predictors they share, and with the same \code{known_sigma}.
}

\value{
An object of class \code{"anova"} (a data frame) with one row per model: the residual degrees of freedom (\code{Res.Df}), the orthogonal objective \eqn{S} (\code{Res.Obj}), the change in degrees of freedom (\code{Df}) and in the objective (\code{Obj}) from the preceding model, and the test statistic with its p-value (\code{F value}, \code{Pr(>F)} or \code{Chisq}, \code{Pr(>Chi)}).
}

\author{
Andrej-Nikolai Spiess
}

\seealso{
\code{\link{onls}}, \code{\link{logLik_o}} for AIC/BIC of non-nested models, \code{\link{profile.onls}} and \code{\link{confint.onls}}, which use the same criterion for single parameters.
}

\examples{
## Enzyme kinetics with errors in the substrate concentration (pipetting, 8\%)
## and in the rate: Michaelis-Menten (h = 1) or cooperative Hill kinetics?
set.seed(42)
S_true <- rep(c(0.25, 0.5, 1, 2, 4, 8, 16, 32), each = 3)
ENZ <- data.frame(S = S_true * exp(rnorm(24, 0, 0.08)),
                  v = 10 * S_true^1.5 / (3^1.5 + S_true^1.5) + rnorm(24, 0, 0.25))

mm    <- onls(v ~ Vm * S / (K + S), data = ENZ, start = list(Vm = 10, K = 3))
hill  <- onls(v ~ Vm * S^h / (K^h + S^h), data = ENZ,
              start = list(Vm = 10, K = 3, h = 1))
hill0 <- onls(v ~ v0 + Vm * S^h / (K^h + S^h), data = ENZ,
              start = list(v0 = 0, Vm = 10, K = 3, h = 1.5))

## Sequential F tests on the orthogonal objective (variances estimated):
## the Hill exponent is needed, a baseline rate is not
anova(mm, hill, hill0)
}

\keyword{models}
\keyword{nonlinear}
\keyword{htest}
