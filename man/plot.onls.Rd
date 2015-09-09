\name{plot.onls}
\alias{plot.onls}

\title{Plotting function for 'onls' objects}

\description{
Plots an orthogonal nonlinear model obtained from \code{\link{onls}}.
}

\usage{
\method{plot}{onls}(x, fitted.nls = TRUE, fitted.onls = TRUE, 
                    segments = TRUE,...)
}

\arguments{
  \item{x}{an object returned from \code{\link{onls}}.}
  \item{fitted.nls}{logical. If \code{TRUE}, the fit from the normal (vertical) nonlinear model is plotted as a blue line for comparison purposes.}
   \item{fitted.onls}{logical. If \code{TRUE}, the fit from the orthogonal nonlinear model is plotted as a red line.}
    \item{segments}{logical. If \code{TRUE}, segments connecting \eqn{(x_i, y_i)} and \eqn{(x_{0i}, y_{0i})} are displayed.}  
  \item{...}{other parameters to \code{\link{plot}}.}       
}

\value{
A plot of the \code{\link{onls}} model.
}
     
\author{
Andrej-Nikolai Spiess
}
        
\examples{
## Quadratic model with 10\% added noise.
## Omitting the "nls" curve,
## display orthogonal segments.
x <- 1:20
y <- 10 + 3*x^2
y <- sapply(y, function(a) rnorm(1, a, 0.1 * a))
DAT <- data.frame(x, y)
mod <- onls(y ~ a + b * x^2, data = DAT, start = list(a = 1, b = 1))
plot(mod, fitted.nls = FALSE)

## Due to different scaling, 
## orthogonality of fit is not evident.
## We need to have equal x/y-scaling for that:
plot(mod, fitted.nls = FALSE, xlim = c(0, 50),
     ylim = c(0, 50))
}

\keyword{models}
\keyword{nonlinear}
