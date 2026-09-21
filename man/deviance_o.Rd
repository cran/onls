\name{deviance_o}
\alias{deviance_o}
\encoding{latin1}

\title{Deviance of the orthogonal residuals}

\description{
Returns the deviance (residual sum of squares) \eqn{\sum_{i=1}^n \hat d_i^2} of the fitted, precision-weighted orthogonal distances from the fitted \code{\link{onls}} model -- the same quantity that \code{\link{onls}} itself minimizes. See 'Details' for what \eqn{\hat d_i} represents in general; a single-predictor, unweighted fit is a special case where this reduces to plain squared Euclidean distance.
}

\usage{
deviance_o(object) 
}

\arguments{
\item{object}{an object returned from \code{\link{onls}}.}
}

\details{
\eqn{\hat d_i} is the fitted, precision-weighted orthogonal distance of observation \eqn{i} (see 'Details' in \code{\link{onls}} for its full definition and the construction of \eqn{Qyy_i}, \eqn{Qx_i} from \code{weights}, \code{sigma_x}, and \code{sigma_y}); it reduces to plain Euclidean distance only for a single-predictor, unweighted fit. \code{deviance_o} returns \eqn{\sum_{i=1}^n \hat d_i^2}, numerically identical to \code{sum(residuals_o(object)^2)} and to the value \code{\link{onls}} itself minimized (\code{object$objective}).
}

\value{
The deviance of the fitted orthogonal (precision-weighted) distances.
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
## See 'onls'.
}   

\keyword{optimize}
\keyword{models}
\keyword{nonlinear}
