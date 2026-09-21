\name{print.onls}
\alias{print.onls}

\title{Printing function for 'onls' objects}

\description{
Provides a printed summary of the converged fit obtained from \code{\link{onls}}.
}

\usage{
\method{print}{onls}(x, ...)
}

\arguments{
  \item{x}{an object returned from \code{\link{onls}}.}
  \item{...}{other parameters for future methods.}       
}

\value{
A printed summary containing the formula, data name, converged parameters, vertical residual sum-of-squares \deqn{\sum_{i=1}^{n} w_i (y_i - \hat{y_i})^2} (prefixed \dQuote{weighted} when \code{weights} are non-constant), orthogonal residual sum-of-squares \deqn{\sum_{i=1}^{n} \hat d_i^2} (see \code{\link{onls}} for the general, precision-weighted definition of \eqn{\hat d_i}; plain Euclidean distance is a special case), a \dQuote{PASSED}/\dQuote{FAILED} message reporting how many fitted points are orthogonal as obtained from \code{\link{check_o}}, and the number of iterations to convergence. 
}
     
\author{
Andrej-Nikolai Spiess
}
        
\examples{
## See 'onls'.
}

\keyword{models}
\keyword{nonlinear}
