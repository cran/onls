\name{deviance_o}
\alias{deviance_o}
\encoding{latin1}

\title{Deviance of the orthogonal residuals}

\description{
Returns the deviance \eqn{\sum_{i=1}^N \| D_i \|^2} of the orthogonal residuals from the fitted \code{\link{onls}} model.
}

\usage{
deviance_o(object) 
}

\arguments{
\item{object}{an object returned from \code{\link{onls}}.}
}

\value{
The deviance of the orthogonal residuals.
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
DNase1 <- subset(DNase, Run == 1)
DNase1$density <- sapply(DNase1$density, function(x) rnorm(1, x, 0.1 * x))
mod <- onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
             data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1))
deviance_o(mod)
}   

\keyword{optimize}
\keyword{models}
\keyword{nonlinear}

