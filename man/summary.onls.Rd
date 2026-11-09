\name{summary.onls}
\alias{summary.onls}

\title{Summary function for 'onls' objects}

\description{
Provides a summary for the parameters of the converged fit, including their standard errors, \eqn{t} values (\eqn{z} values for known variances) and p values, together with residual standard errors with respect to both the vertical and the orthogonal residuals. See 'Details' for how parameters held \code{fixed} in the original \code{\link{onls}} call, and \code{known_sigma}-based scaling, are handled.
}

\usage{
\method{summary}{onls}(object, correlation = FALSE, symbolic.cor = FALSE, ...)
}

\arguments{
  \item{object}{an object returned from \code{\link{onls}}.}
  \item{correlation}{logical. If \code{TRUE}, the correlation matrix of the estimated (non-\code{fixed}) parameters is returned and printed.}
  \item{symbolic.cor}{logical. If \code{TRUE}, print the correlations in a symbolic form.}
  \item{...}{further arguments passed to or from other methods.}       
}

\details{
\bold{Fixed parameters.} If \code{object} was fitted with some parameters held \code{fixed} (see \code{\link{onls}}), those parameters were never estimated: their \code{Estimate} is still shown (the fixed value used), but \code{Std. Error}, \code{t value} and \code{Pr(>|t|)} are reported as \code{NA} rather than \code{0}/\code{Inf}/a spuriously small p-value, matching R's own convention for non-estimated coefficients (e.g. aliased terms in \code{\link{summary.lm}}). The residual degrees of freedom in \code{df} count only the \emph{free} (non-fixed) parameters, so that \code{df[1] + df[2]} equals the number of observations with non-zero weight, as in \code{\link{summary.lm}}/\code{\link{summary.nls}}. If \code{correlation = TRUE}, fixed parameters are dropped from the returned correlation matrix entirely (rather than padded with \code{NA} rows/columns), since correlation with a held-fixed constant is not a meaningful quantity.

\bold{Known variances.} If the variances were treated as known (\code{known_sigma = TRUE} in \code{\link{onls}}), the standard errors are not estimated from the residuals, and the test statistics are compared with the normal distribution (\code{z value}, \code{Pr(>|z|)}); otherwise with the \eqn{t} distribution with \code{df_resid} degrees of freedom. \code{\link{confint.onls}} and \code{\link{predict.onls}} use the same rule.

\bold{cov.unscaled.} \code{object$vcov} was already scaled by \code{\link{onls}} itself: not at all if \code{known_sigma = TRUE}, or by the \emph{orthogonal}-distance reduced chi-square (\code{object$reduced_chisq}) if \code{known_sigma = FALSE}. \code{cov.unscaled} undoes exactly that scaling (dividing by \code{object$reduced_chisq} only in the latter case), so that \code{vcov(object) == cov.unscaled * sigma^2} holds with the appropriate \code{sigma}, matching the \code{\link{summary.lm}}/\code{\link{summary.nls}} convention. This is \emph{not} the same quantity as the \emph{vertical}-residual variance (\code{sigmaONLS^2}): the vertical and orthogonal residual sums of squares generally differ, even for an unweighted fit, so \code{cov.unscaled} must not be confused with, or derived from, \code{sigmaONLS}.
}

\value{
A list of class \code{"summary.onls"}, with components:
\item{formula}{the model formula.}
\item{residONLS}{the vertical residuals (observed minus fitted response at the observed predictors), as \code{\link{residuals}}.}
\item{sigmaONLS}{residual standard error of the \emph{vertical} distances.}
\item{sigma_o}{residual standard error of the \emph{orthogonal} (precision-weighted) distances.}
\item{df}{a length-2 vector, the number of free (non-\code{fixed}) parameters and the residual degrees of freedom; \code{sum(df)} equals the number of observations with non-zero weight.}
\item{cov.unscaled}{the unscaled covariance matrix of the free parameters; see 'Details'.}
\item{call}{the matched call.}
\item{convInfo}{convergence information from \code{\link{onls}}.}
\item{control}{the control settings used, see \code{\link{onls}}.}
\item{na.action}{information on the handling of \code{NA}s.}
\item{coefficients, parameters}{(identical) matrices of \code{Estimate}, \code{Std. Error}, \code{t value} and \code{Pr(>|t|)} (\code{z value} and \code{Pr(>|z|)} if the variances are known), one row per parameter (including \code{fixed} ones, with \code{NA} in the last three columns for those); see 'Details'.}
\item{known_sigma}{as in \code{object}.}
\item{reduced_chisq}{as in \code{object}.}
\item{correlation}{(only if \code{correlation = TRUE}) the correlation matrix of the \emph{free} parameters only; see 'Details'.}
\item{symbolic.cor}{(only if \code{correlation = TRUE}) as supplied.}
}
     
\author{
Andrej-Nikolai Spiess
}
        
\examples{
## See 'onls'.
}

\keyword{models}
\keyword{nonlinear}
