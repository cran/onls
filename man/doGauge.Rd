\name{doGauge}
\alias{doGauge}

\title{Gauge onls against published, closed-form, ODRPACK and NIST reference values}

\description{
Runs the demo \code{gauge} of the package, which fits the reference problems of orthogonal distance regression and compares every result of \code{\link{onls}} (parameters, standard errors, objective) with a reference value, and shows or writes its complete output (code and results). The result of each comparison (\code{PASS}, \code{NOTE} or \code{DIFF}) is highlighted. A summary of all comparisons is printed and returned.
}

\usage{
doGauge(output = c("console", "html", "pdf"), outdir = tempdir(), fontsize = 14,
        silent = TRUE, file = system.file("demo", "gauge.R", package = "onls"))
}

\arguments{
  \item{output}{\code{"console"} (default): the output is printed; \code{"html"} or \code{"pdf"}: it is written to the file \code{onls-gauge.html} or \code{onls-gauge.pdf} in \code{outdir}.}
  \item{outdir}{the directory in which the HTML or PDF file is written. The default is a temporary directory; use, \emph{e.g.}, \code{getwd()} to keep the file.}
  \item{fontsize}{the font size of the output in pixels (HTML); for PDF, the size is scaled accordingly and reduced if necessary, so that the longest line fits the page width.}
  \item{silent}{logical. If \code{TRUE}, messages and warnings issued while the demo runs are suppressed.}
  \item{file}{the demo file; only needed to run a modified copy of the demo.}
}

\details{
The demo contains five kinds of comparison (column 2 of the output):
\itemize{
  \item \code{published}: values printed in the cited publications: the ODRPACK User's Reference Guide (example 2.C.i), Algorithm 676 (ODRPACK), the bounded example of ODRPACK95 (TOMS 869), Pearson's data with York's weights (York 1966; values of the York et al. 2004 equations, as reported by the R package 'bfsl'), Daeron & Vermeesch (2024, Table 3) and the XLSTAT Deming regression example. The tolerance is half a unit of the last printed decimal; for results of iterative ODRPACK runs, the convergence of that run is added (\eqn{10^{-3}} of the standard error for parameters, a relative \eqn{10^{-4}} for standard errors and \eqn{10^{-7}} for the objective).
  \item \code{closed form}: exact solutions computed in the demo: total least squares and generalized (correlated-error) total least squares by singular value decomposition (Golub & Van Loan 1980), and Deming regression; relative tolerance \eqn{10^{-6}}.
  \item \code{ODRPACK}: ODRPACK 2.01 (Boggs et al. 1989) through \code{scipy.odr} 1.18.1, with the same data, model, weights and starting values, run to convergence (\code{sstol = partol = 1e-14}), for problems without published orthogonal-regression results: the example with separate response and predictor weights (ODRPACK's \code{WE}, \code{WD}), and the NIST problems Chwirut2, Gauss1, ENSO and the Rat43 data (Gompertz model) as orthogonal regressions with unit precisions. Parameters are compared to \eqn{10^{-3}} of their standard error, standard errors to a relative \eqn{10^{-4}}, the objective to a relative \eqn{10^{-7}}.
  \item \code{consistency}: properties that must hold (the functions around the fit, and the agreement of the two algorithms); \code{PASS} if the property holds.
  \item \code{NIST}: the certified values of the NIST StRD nonlinear regression data sets Chwirut2, Gauss1, ENSO and Rat43, which refer to ordinary least squares. With a negligible \code{sigma_x}, the orthogonal fit becomes ordinary least squares, so \code{onls} must reproduce the certified parameters (relative tolerance \eqn{10^{-6}}) and residual sums of squares (\eqn{10^{-9}}).
}
Results:
\itemize{
  \item \code{PASS}: agreement within the tolerance.
  \item \code{NOTE}: a documented discrepancy of the reference, explained in a comment of the demo (the ODRPACK95 reference run reports a higher objective than the bounded optimum; a truncated instead of rounded published value).
  \item \code{DIFF}: a deviation that is neither within the tolerance nor documented. The expected result is \code{DIFF = 0}.
}
The demo can also be run directly with \code{demo("gauge", package = "onls")}. It replaces the reference examples of \code{\link{onls}} and the function \code{NIST} of earlier versions.
}

\value{
Invisibly, a list with
\item{file}{the path of the file written (\code{NULL} for \code{output = "console"}).}
\item{results}{a data frame with one row per comparison: \code{item}, \code{kind} (of reference), \code{reference}, \code{onls}, \code{diff}, \code{tol} and \code{result}.}
\item{summary}{a table with the numbers of \code{PASS}, \code{NOTE} and \code{DIFF}.}
}

\author{
Andrej-Nikolai Spiess
}

\references{
ALGORITHM 676 ODRPACK: Software for Weighted Orthogonal Distance Regression.\cr
Boggs PT, Donaldson JR, Byrd RH and Schnabel RB.\cr
ACM Trans Math Soft (1989), \bold{15}, 348-364.\cr
\doi{10.1145/76909.76913}

User's Reference Guide for ODRPACK Version 2.01. Software for Weighted Orthogonal Distance Regression.\cr
Boggs PT, Byrd RH, Rogers JE and Schnabel RB.\cr
NISTIR (1992), \bold{4834}: 1-113.

Algorithm 869: ODRPACK95: A weighted orthogonal distance regression code with bound constraints.\cr
Zwolak JW, Boggs PT and Watson LT.\cr
ACM Trans Math Soft (2007), \bold{33}: 27.\cr
\doi{10.1145/1268776.1268782}

Least squares fitting of a straight line.\cr
York D.\cr
Canadian Journal of Physics (1966), \bold{44}: 1079-1086.\cr
\doi{10.1139/p66-090}

Unified equations for the slope, intercept, and standard errors of the best straight line.\cr
York D, Evensen NM, Martinez ML and De Basabe Delgado J.\cr
American Journal of Physics (2004), \bold{72}: 367-375.

Omnivariant generalized least squares regression.\cr
Daeron M and Vermeesch P.\cr
Chemical Geology (2024), \bold{647}: 121881.\cr
\doi{10.1016/j.chemgeo.2023.121881}

An analysis of the Total Least Squares problem.\cr
Golub GH and Van Loan CF.\cr
SIAM Journal on Numerical Analysis (1980), \bold{17}: 883-893.\cr
\doi{10.1137/0717073}

Statistical Reference Datasets (StRD), nonlinear regression.\cr
National Institute of Standards and Technology.\cr
\url{https://www.itl.nist.gov/div898/strd/nls/nls_main.shtml}
}

\seealso{
\code{\link{onls}}, \code{\link{check_o}}.
}

\examples{
\donttest{
res <- doGauge()
res$summary
## HTML file in the working directory:
## doGauge(output = "html", outdir = getwd())
}
}

\keyword{models}
\keyword{nonlinear}
