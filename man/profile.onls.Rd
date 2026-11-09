\name{profile.onls}
\alias{profile.onls}
\alias{plot.profile.onls}
\alias{confint.profile.onls}

\title{Profiles of the orthogonal-distance objective}

\description{
Profiles the minimized orthogonal-distance objective of an \code{\link{onls}} fit for each free parameter, and plots the profiles or turns them into confidence intervals.
}

\usage{
\method{profile}{onls}(fitted, which = NULL, maxpts = 100, alphamax = 0.01,
        delta.t = cutoff/5, ...)
\method{plot}{profile.onls}(x, absVal = TRUE, level = 0.95, ...)
\method{confint}{profile.onls}(object, parm, level = 0.95, ...)
}

\arguments{
  \item{fitted}{an object returned from \code{\link{onls}}.}
  \item{which}{parameters to profile, as names or indices; default all free parameters. Fixed parameters are not profiled.}
  \item{maxpts}{maximum number of points in each direction of a profile.}
  \item{alphamax}{the profile is extended in each direction until \eqn{|\tau|} exceeds the \eqn{1 - } \code{alphamax} quantile (see 'Details').}
  \item{delta.t}{step size of the profile on the \eqn{\tau} scale; \code{cutoff} is that quantile.}
  \item{x, object}{an object returned from \code{profile.onls}.}
  \item{absVal}{logical. If \code{TRUE}, \eqn{|\tau|} is plotted, otherwise \eqn{\tau}.}
  \item{level}{confidence level for the reference lines of the plot and for the intervals.}
  \item{parm}{parameters for which intervals are computed, as names or indices.}
  \item{...}{further arguments, passed to \code{\link{plot}} by the plot method and ignored otherwise.}
}

\details{
For parameter \eqn{\theta_j}, the parameter is held fixed at values on both sides of its estimate, and all other free parameters and all foot points are re-optimized by \code{\link{onls}} (using the data, weights and settings stored in \code{fitted}, so that the profile does not depend on objects in the calling environment). With \eqn{\hat S} the minimized objective and \eqn{S(\theta_j)} the constrained minimum, the profile statistic is
\deqn{\tau(\theta_j) = \mathrm{sign}(\theta_j - \hat\theta_j) \sqrt{\left(S(\theta_j) - \hat S\right)/s^2}}
with \eqn{s^2 = \hat S/\mathrm{df}} when the variances are estimated and \eqn{s^2 = 1} when they are known (\code{known_sigma = TRUE}). The grid is extended until \eqn{|\tau|} exceeds \eqn{\sqrt{F_{1-\alpha}(1, \mathrm{df})}} (estimated variances) or \eqn{\sqrt{\chi^2_{1-\alpha}(1)}} (known variances), with \eqn{\alpha} = \code{alphamax}, or a failed refit, a flat profile, a \code{lower}/\code{upper} bound (at which the profile is evaluated last) or \code{maxpts} stops it; the reason is stored in the \code{stop} element of the \code{"parameters"} attribute of each profile. The stepping follows \code{\link[stats]{profile.nls}}. Models with one free parameter are profiled by refits in which only the foot points are optimized. With estimated variances, an exact fit (no residual degrees of freedom) cannot be profiled and gives an error.

\code{confint} uses the part of the profile on which \eqn{\tau} increases steadily away from the estimate, interpolates it (a monotone spline of \eqn{\tau} over the parameter, inverted) at \eqn{\tau = \pm q}, with \eqn{q} the normal (known variances) or \eqn{t} quantile; \code{\link{confint.onls}} with \code{method = "profile"} calls it.
}

\value{
\code{profile.onls}: a list with one data frame per profiled parameter, with components \code{tau} and \code{par.vals} (the parameter values along the profile), of class \code{c("profile.onls", "profile")}. The attributes \code{"original.fit"} and \code{"summary"} hold the fit and its summary.

\code{confint.profile.onls}: a matrix of lower and upper limits. \code{plot.profile.onls}: \code{NULL}, invisibly.
}

\author{
Andrej-Nikolai Spiess
}

\seealso{
\code{\link{confint.onls}}.
}

\examples{
set.seed(123)
DNase1 <- subset(DNase, Run == 1)
DNase1$density <- sapply(DNase1$density, function(x) rnorm(1, x, 0.1 * x))
mod1 <- onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
             data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1))
prof <- profile(mod1)
plot(prof)
confint(prof)
}

\keyword{models}
\keyword{nonlinear}
