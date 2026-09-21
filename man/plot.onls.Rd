\name{plot.onls}
\alias{plot.onls}

\title{Plotting function for 'onls' objects}

\description{
Plots orthogonal nonlinear models obtained from \code{\link{onls}}: the observed data, the fitted \code{\link{onls}} curve (and, optionally, the ordinary \code{\link{nls}} warm-start curve), and segments connecting each observation \eqn{(x_i, y_i)} to its foot point \eqn{(x_{0i}, y_{0i})}. 
}

\usage{
\method{plot}{onls}(x, panel = NULL, fitted.nls = TRUE, fitted.onls = TRUE, asp = TRUE, 
                    segments = TRUE, npoints = 200, nmesh = 25, ...)
}

\arguments{
\item{x}{an object returned from \code{\link{onls}}.}
\item{panel}{for multivariate models only (more than one predictor): an integer index or predictor name selecting a single predictor to plot as one full-size plot, with the same treatment (including \code{asp = 1} by default) as a single-predictor model. If \code{NULL} (default), all predictors are plotted together as a grid of smaller panels instead (see 'Details'). Ignored, with a warning, for single-predictor models.}
\item{fitted.nls}{logical. If \code{TRUE}, the fit from the normal (vertical) nonlinear model is plotted as a blue line for comparison purposes.}
\item{fitted.onls}{logical. If \code{TRUE}, the fit from the orthogonal nonlinear model is plotted as a red line.}
\item{asp}{logical. If \code{TRUE}, plots are generated in a 1:1 axis ratio for exact display of the orthogonality.}  
\item{segments}{logical. If \code{TRUE}, segments connecting \eqn{(x_i, y_i)} and \eqn{(x_{0i}, y_{0i})} are displayed.}  
\item{npoints}{number of points used to draw the smooth \code{fitted.onls}/\code{fitted.nls} curve(s).}
\item{nmesh}{number of \code{nmesh x nmesh} divisions used to draw the 3D mesh when two predictors are used.}
\item{...}{other parameters to \code{\link{plot}}, such as \code{xlim}, \code{ylim} or \code{asp}.}       
}

\details{
There are three plot types, depending on the model and on \code{panel}:

\bold{Single-predictor plot} => \emph{p} = 1; single-predictor models, or a multivariate model with \code{panel} set to one predictor):\cr
a single, full-size plot of \eqn{y} against that predictor, with \code{asp = 1} used by default (unless overridden via \code{...}). 
Under \code{asp = 1}, one data unit in \eqn{x} is rendered the same physical length as one data unit in \eqn{y}, so for an \emph{unweighted} \code{\link{onls}} fit 
the segments connecting \eqn{(x_i,y_i)} to \eqn{(x_{0i},y_{0i})} genuinely appear at right angles to the fitted curve -- a visual check of what \code{\link{check_o}} verifies numerically. 
If only one of \code{xlim}/\code{ylim} is supplied via \code{...}, the other is automatically chosen (via \code{\link{extendrange}} on the data that falls inside the supplied range) so that both remain sensibly matched to the data before \code{asp = 1} does its own adjustment; 
supplying only one of the two is therefore usually enough to zoom in on a region of interest.

\bold{3D plot} => \emph{p} = 2; bivariate models:\cr
An rgl-based 3D plot with x1/x2/z setup, and where the segments connect the estmated values to the foot points on the 3D surface.
The advantage of this plot is the rotation- and zoom-ability to inspect orthogonality of the points.

\bold{Multivariate grid} => \emph{p > 2}; multivariate models:\cr
One partial-dependence panel per predictor, each showing \eqn{y} against one predictor with the \emph{other} predictors held fixed at their mean foot-point value. 
\code{asp} is \emph{not} forced to \code{1} here (unless explicitly supplied via \code{...}): panels are shrunk by the \code{mfrow} layout, and forcing \code{asp = 1} 
in a small panel can require aggressively re-expanding the axes to preserve a strict 1:1 unit scale, which would silently override any \code{xlim}/\code{ylim} requested for that panel. 
These grid panels are also not a literally faithful orthogonality check to begin with: since the curve fixes the \emph{other} predictors at a shared mean value while a given point's own segment 
endpoint uses its own foot point in every dimension, the two generally do not coincide exactly on the drawn curve, so segments should not be expected to look exactly perpendicular even for an unweighted fit. 
Use \code{\link{check_o}} for a reliable, per-observation, per-axis orthogonality check rather than eyeballing these panels; use \code{panel} to get a single full-size, \code{asp = 1} plot for one predictor at a time instead.
}

\value{
A plot of the \code{\link{onls}} model (single-predictor plot, single chosen \code{panel}, multivariate grid or 3D).
}
     
\author{
Andrej-Nikolai Spiess
}
        
\examples{
## 1a. Quadratic model with 10\% added noise.
set.seed(123)
x <- 1:20
y <- 10 + 3*x^2
y <- y + rnorm(20, 0, 50)
DAT <- data.frame(x, y)
mod1 <- onls(y ~ a + b * x^2, data = DAT, start = list(a = 10, b = 3), extend = c(0.2, 0))
plot(mod1)

## 1b. Zooming in on a region: supplying xlim alone is enough --
## a matching ylim is chosen automatically before asp = 1 is applied.
## Need to set asp = FALSE, but a 1-to-1 aspect is not exactly possible.
plot(mod1, fitted.nls = FALSE, xlim = c(0, 10), asp = FALSE)

## 2. Half-dome 3D example with rgl plot
set.seed(123)
n <- 60
r_true <- 6
ang <- runif(n, 0, 2 * pi)
rad <- sqrt(runif(n, 0, 0.55)) * r_true
x1  <- rad * cos(ang)
x2  <- rad * sin(ang)
z <- sqrt(r_true^2 - x1^2 - x2^2) + rnorm(n, 0, 0.15)
x1 <- x1 + rnorm(n, 0, 0.1)
x2 <- x2 + rnorm(n, 0, 0.1)
DAT <- data.frame(x1 = x1, x2 = x2, z = z)
maxrad <- max(sqrt(x1^2 + x2^2))
mod2 <- onls(z ~ sqrt(r^2 - x1^2 - x2^2), data = DAT,
                 start = list(r = r_true),
                 sigma_x = c(0.1, 0.1), sigma_y = 0.15,
                 lower = maxrad * 1.05, upper = 100)
check_o(mod2) # all orthogonal to dome surface
plot(mod2)    # should render as a visibly ROUND dome

\donttest{
## 3. Multivariate setup with 4 predictors
set.seed(123)
n  <- 60
x1 <- runif(n, 0, 10)
x2 <- runif(n, 0, 5)
x3 <- runif(n, 2, 10)   
x4 <- runif(n, 2, 10)  

b0 <- 2; b1 <- 0.8; b2 <- 0.5; b3 <- 2; b4 <- 3
z_true <- b0 + b1 * x1 + b2 * x2^2 + b3 * sqrt(x3) + b4 * log(x4 + 1)
z <- z_true + rnorm(n, 0, 0.5)

sd_x <- c(0.3, 0.2, 0.3, 0.3)
x1 <- x1 + rnorm(n, 0, sd_x[1])
x2 <- x2 + rnorm(n, 0, sd_x[2])
x3 <- x3 + rnorm(n, 0, sd_x[3])
x4 <- x4 + rnorm(n, 0, sd_x[4])

DAT <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, z = z)
mod3 <- onls(z ~ b0 + b1 * x1 + b2 * x2^2 + b3 * sqrt(x3) + b4 * log(x4 + 1),
            data = DAT, start = list(b0 = 1, b1 = 1, b2 = 1, b3 = 1, b4 = 1),
            sigma_x = sd_x, sigma_y = 0.5)

summary(mod3)      
check_o(mod3)       
plot(mod3)
}
}

\keyword{models}
\keyword{nonlinear}
