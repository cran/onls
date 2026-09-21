\name{residuals_o}
\alias{residuals_o}
\encoding{latin1}

\title{The orthogonal residuals}

\description{
Returns a vector of the fitted, precision-weighted orthogonal distances \eqn{\hat d_i} from the fitted \code{\link{onls}} model, in the original observation order -- the same per-observation quantities that are squared and summed by \code{\link{deviance_o}}. See 'Details' for what \eqn{\hat d_i} represents in general; a single-predictor, unweighted fit is a special case where this reduces to plain Euclidean distance.
}

\usage{
residuals_o(object) 
}

\arguments{
\item{object}{an object returned from \code{\link{onls}}.}
}

\details{
\eqn{\hat d_i} is the fitted, precision-weighted orthogonal distance of observation \eqn{i} (see 'Details' in \code{\link{onls}} for its full definition and the construction of \eqn{Qyy_i}, \eqn{Qx_i} from \code{weights}, \code{sigma_x}, and \code{sigma_y}); it reduces to plain Euclidean distance only for a single-predictor, unweighted fit. Values are returned in the original observation order regardless of the number of predictors, and any rows omitted via \code{na.action} are correctly reinserted as \code{NA} (via \code{\link{napredict}}) at their original positions.
}

\value{
A vector of the fitted orthogonal (precision-weighted) distances.
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
