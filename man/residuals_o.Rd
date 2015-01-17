\name{residuals_o}
\alias{residuals_o}
\encoding{latin1}

\title{The orthogonal residuals}

\description{
Returns a vector with the orthogonal residuals \eqn{\| D_i \|} from the fitted \code{\link{onls}} model.
}

\usage{
residuals_o(object) 
}

\arguments{
\item{object}{an object returned from \code{\link{onls}}.}
}

\value{
The orthogonal residuals.
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
DNase1 <- subset(DNase, Run == 1)
DNase1$density <- sapply(DNase1$density, function(x) rnorm(1, x, 0.1 * x))
mod <- onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
             data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1))
residuals_o(mod)
}   

\keyword{optimize}
\keyword{models}
\keyword{nonlinear}

