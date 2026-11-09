\name{confint.onls}
\alias{confint.onls}

\title{Confidence intervals for 'onls' model parameters}

\description{
Computes confidence intervals for the free parameters of an \code{\link{onls}} model, from the profile of the orthogonal-distance objective (default), from the standard errors (Wald), or by a nonparametric bootstrap.
}

\usage{
\method{confint}{onls}(object, parm, level = 0.95, method = c("profile", "wald", "boot"),
        k = 200, max_tries = 10 * k, verbose = FALSE, ...)
}

\arguments{
  \item{object}{an object returned from \code{\link{onls}}.}
  \item{parm}{parameters for which intervals are computed, as names or indices; default all free (non-\code{fixed}) parameters. Fixed parameters are dropped with a warning.}
  \item{level}{the confidence level required.}
  \item{method}{\code{"profile"} (default), \code{"wald"} or \code{"boot"}; see 'Details'.}
  \item{k}{for \code{method = "boot"}: the number of successful bootstrap fits.}
  \item{max_tries}{for \code{method = "boot"}: the maximum number of refits attempted; if fewer than \code{k} succeed, a warning gives the count.}
  \item{verbose}{for \code{method = "boot"}: logical, print a progress line.}
  \item{...}{for \code{method = "profile"}: further arguments passed to \code{\link{profile.onls}}, such as \code{maxpts}.}
}

\details{
The profile and Wald intervals use the quantile \eqn{q} of the normal distribution when the variances are known (\code{known_sigma = TRUE} in \code{\link{onls}}) and of the \eqn{t} distribution with \code{df_resid} degrees of freedom otherwise.

\bold{\code{"profile"}.} Each parameter is held fixed on a grid around its estimate while all other free parameters and all foot points are re-optimized (\code{\link{profile.onls}}). With \eqn{S} the minimized orthogonal objective, the profile statistic is \eqn{\tau = \mathrm{sign}(\theta_j - \hat\theta_j)\sqrt{(S(\theta_j) - \hat S)/s^2}}, with \eqn{s^2 = \hat S/\mathrm{df}} (estimated variances) or \eqn{s^2 = 1} (known variances). The limits are where \eqn{\tau = \pm q}, found by interpolating the profile. Profile intervals follow the shape of the objective and can be asymmetric. Only the part of the profile on which \eqn{\tau} increases steadily away from the estimate is used. A limit that the profile does not reach is returned as \code{NA}, with a warning that says why: most often the profile flattens out because the data do not bound the parameter on that side, or a \code{lower}/\code{upper} bound lies inside the interval.

\bold{\code{"wald"}.} \eqn{\hat\theta_j \pm q \cdot SE_j}, with the standard errors of \code{\link{vcov}}. Instant and accurate for well-determined parameters, but symmetric by construction, and too narrow when the objective is strongly asymmetric in a parameter.

\bold{\code{"boot"}.} Nonparametric case bootstrap. In each replicate, the rows of the data used in the fit (after \code{subset} and \code{na.action}) are sampled with replacement, together with their row-aligned \code{weights}, \code{sigma_x} and \code{sigma_y}, and the model is refitted with all other settings of the original call, starting from the original estimates. Fits that do not converge are discarded; if all points of the original fit were orthogonal, refits with non-orthogonal points are discarded as well. The limits are the empirical quantiles of the \code{k} successful fits. At most \code{max_tries} refits are attempted. Case resampling is unreliable for very small \eqn{n}.
}

\value{
A matrix with one row per parameter and columns giving the lower and upper confidence limits. For \code{method = "boot"}, the attribute \code{"n_boot"} holds the number of successful bootstrap fits.
}
     
\author{
Andrej-Nikolai Spiess
}

\seealso{
\code{\link{profile.onls}}, \code{\link{vcov}}.
}
        
\examples{
set.seed(123)
DNase1 <- subset(DNase, Run == 1)
DNase1$density <- sapply(DNase1$density, function(x) rnorm(1, x, 0.1 * x))
mod1 <- onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
             data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1))
confint(mod1)                      # profile
confint(mod1, method = "wald")     # symmetric, too narrow for Asym
\donttest{
confint(mod1, method = "boot", k = 50) # just time optimization, consider k >= 200.
}
}

\keyword{models}
\keyword{nonlinear}
