\name{onls}
\alias{onls}
\encoding{UTF-8}

\title{Weighted Orthogonal Nonlinear Least-Squares Regression (Orthogonal Distance Regression)}

\description{
  Fits nonlinear orthogonal-distance regression (ODR, also called errors-in-variables) models in the manner of ODRPACK (Boggs, Byrd and Schnabel, 1987; Boggs et al., 1992): the model parameters and one foot-point correction per observation are estimated \emph{simultaneously} by minimizing the same weighted sum of squared response and predictor residuals that ODRPACK minimizes, using a Levenberg-Marquardt algorithm on the joint problem. 
  In contrast to ordinary nonlinear least squares, \code{onls()} allows measurement error in both predictors and response variables. 
  Predictor uncertainty may be specified through predictor-specific standard deviations, observation-specific error structures, or full covariance matrices allowing correlated predictor errors. Response uncertainty may be supplied globally or on a per-observation basis and can be combined with observation weights. 
  The algorithm supports single- and multi-predictor nonlinear models, parameter bounds, fixed parameters, and covariance estimation for the fitted parameters that is equivalent to the ODRPACK (Gauss-Newton) covariance.
}

\usage{
  onls(formula, data, start = NULL, weights = NULL, sigma_x = NULL, sigma_y = 1, 
       known_sigma = NULL, extend = NULL, window = NULL, control = list(), 
       lower = NULL, upper = NULL, fixed = NULL, subset = NULL, na.action = NULL, 
       trace = FALSE, algorithm = c("auto", "structured", "dense"))
}

\arguments{
  \item{formula}{A two-sided nonlinear model \code{\link{formula}} describing the response as a nonlinear function of one or more predictors and unknown model parameters. The response must be a plain variable name (no transformation of the left-hand side).}
  \item{data}{A data frame containing the variables appearing in \code{formula}.}
  \item{start}{A named list (or named vector) of starting values for all model parameters.}
  \item{weights}{Optional non-negative observation weights. As in \code{\link{nls}}, \code{weights}, \code{subset}, \code{sigma_x} and \code{sigma_y} are evaluated in \code{data} first and then in the calling environment, so that columns of \code{data} can be named directly (\emph{e.g.} \code{weights = w}, \code{subset = conc > 0.1}; \code{NA} in a logical \code{subset} drops the observation). The weights are incorporated into the orthogonal-distance criterion through the response precision matrix and therefore affect both the foot-point corrections and the parameter estimates. When supplied, \code{weights} must have length equal to the number of observations.}
  \item{sigma_x}{Specification of predictor measurement error, as a standard deviation (not a variance or a weight). Accepts \code{NULL}, a scalar, a length-\eqn{p} vector, a \eqn{p \times p} covariance matrix, or an \eqn{n \times p} matrix of observation-specific standard deviations (a length-\eqn{n} vector is accepted as shorthand for the latter when \eqn{p = 1}). See 'Details'.}
  \item{sigma_y}{Response measurement error, as a standard deviation. May be either a single positive value applied to all observations or a length-\eqn{n} vector of observation-specific standard deviations.}
  \item{known_sigma}{Logical indicating whether the supplied (or default) \code{sigma_x}/\code{sigma_y} should be treated as known measurement-error standard deviations rather than as estimated from the data. If \code{NULL} (default), it is set to \code{TRUE} if \code{sigma_x} is supplied or \code{sigma_y} is explicitly passed by the user (regardless of its value), and \code{FALSE} otherwise. Standard errors are scaled by the reduced chi-square statistic (see \code{reduced_chisq} in 'Value') only when \code{known_sigma = FALSE}; when \code{known_sigma = TRUE} the supplied variances are assumed to already be on the correct absolute scale and no such scaling is applied. If \code{known_sigma} is \code{TRUE} but \code{sigma_x} is not given, a warning points out that the predictor errors are then taken to have standard deviation 1 in the units of each predictor.}
  \item{extend}{Optional numeric vector of length one or two. If \code{NULL} (default), the foot-point corrections are unbounded, as in ODRPACK. If supplied (single-predictor models only), the foot points are confined to a range extended beyond the observed predictor range by the given fractions (lower, upper; a single value is used for both ends). See 'Details'.}
  \item{window}{Optional integer window width, used only if \code{extend} is supplied, for single-predictor models with \eqn{n > 25}: restricts the foot point of each observation to the neighborhood of the corresponding sorted observations. See 'Details'.}
  \item{control}{Optional list of optimization settings. Recognized elements are \code{outer_max}, the \emph{total} budget of Levenberg-Marquardt iterations (default \code{5000}; with \code{algorithm = "dense"} spent in restarts of at most 1024 iterations each, see 'Details'); \code{ftol} and \code{ptol}, the relative-reduction and relative-step convergence tolerances (as in \code{\link[minpack.lm]{nls.lm.control}}; default \code{1e-10} each); \code{diag}, an optional length-\eqn{q} vector of positive parameter scaling factors (foot-point corrections are given the scale 1; if omitted, MINPACK's adaptive scaling is used); and \code{lower}/\code{upper}, alternatives to the arguments of the same name.}
  \item{lower}{Optional vector of lower bounds for the model parameters, of length \eqn{q} (same ordering as \code{start}, fixed parameters included).}
  \item{upper}{Optional vector of upper bounds for the model parameters, of length \eqn{q} (same ordering as \code{start}, fixed parameters included).}
  \item{fixed}{Optional logical vector with the same length and ordering as \code{start}. Parameters marked \code{TRUE} are held fixed at their supplied starting values throughout optimization. If all parameters are fixed, only the foot points are optimized (as used by \code{\link{profile.onls}} for one-parameter models).}
  \item{subset}{Optional logical or index expression, evaluated in \code{data}, selecting the observations to be used for fitting. Row-aligned \code{weights}, \code{sigma_x} and \code{sigma_y} are subset accordingly.}
  \item{na.action}{Function indicating how missing values should be handled. Row-aligned \code{weights}, \code{sigma_x} and \code{sigma_y} follow the removal of observations.}
  \item{trace}{Logical. If \code{TRUE}, progress messages from the NLS warm start and the orthogonal-distance optimization (iterations or restarts) are printed.}
  \item{algorithm}{The solver for the joint problem; both minimize the same objective. \code{"auto"} (default) uses \code{"structured"} for more than 100 observations or when parameter bounds (\code{lower}/\code{upper}) are given, and \code{"dense"} otherwise. \code{"structured"} uses the block structure of the problem, as ODRPACK does, and its cost per iteration grows linearly with \eqn{n}. \code{"dense"} passes the full joint problem to \code{\link[minpack.lm]{nls.lm}} (the solver of \pkg{onls} 0.2), with a cost per iteration growing with \eqn{n^3}. Fits with the foot-point bounds \code{extend}/\code{window} always use \code{"dense"}. See 'Details'.}
}

\details{
\bold{Model and objective.} Assume a nonlinear model \eqn{y = f(x, \theta)}, where \eqn{x} is a vector of \eqn{p} predictors and \eqn{\theta} is a vector of \eqn{q} unknown model parameters.
Unlike ordinary nonlinear least squares, which assumes error only in the response, ODR assumes that both the response and the predictors are observed with error. 
For observation \eqn{i}, let \eqn{x_i} be the observed predictor vector, \eqn{y_i} the observed response, \eqn{\delta_i} the unknown correction to the predictors (so that \eqn{\xi_i = x_i + \delta_i} is the \emph{foot point} on the model surface),
\eqn{Qyy_i} the response precision and \eqn{Qx_i} the predictor precision matrix (both constructed below). 
The parameters and all \eqn{n} corrections are estimated jointly by minimizing
\deqn{S(\theta, \delta_1, \dots, \delta_n) = \sum_{i=1}^{n} \left[ Qyy_i \left(y_i - f(x_i + \delta_i, \theta)\right)^2 + \delta_i^T Qx_i \, \delta_i \right]}
This is the explicit ODR problem of ODRPACK, with \eqn{Qyy_i} playing the role of the response weights (\code{WE}) and \eqn{Qx_i} the role of the predictor weights (\code{WD}). It has \eqn{q_{free} + np} unknowns, where \eqn{q_{free}} is the number of free (non-\code{fixed}) parameters.
For fixed \eqn{\theta}, minimizing \eqn{S} over \eqn{\delta_i} alone gives the weighted squared orthogonal distance of observation \eqn{i} to the model surface,
\deqn{d_i^2 = Qyy_i \left[y_i - f(\hat\xi_i, \theta) \right]^2 + (\hat\xi_i - x_i)^T Qx_i (\hat\xi_i - x_i)}
so that the estimate \eqn{\hat\theta} is the minimizer of \eqn{\sum_i d_i^2}; the \eqn{d_i} are returned as \code{orth_dist} (\code{\link{residuals_o}}), and \eqn{\sum_i d_i^2} as \code{objective}.
  
\bold{Predictor precision (\code{sigma_x}):}\cr
The predictor precision \eqn{Qx_i} used above is constructed from \code{sigma_x}, which is always supplied on the standard-deviation scale:\cr
\code{NULL}: unit variance for every predictor, \eqn{Qx_i = I_p} for all \eqn{i}\cr
scalar \eqn{\sigma}: isotropic, \eqn{Qx_i = I_p/\sigma^2} for all \eqn{i}\cr
length-\eqn{p} vector \eqn{(\sigma_{x1}, \dots, \sigma_{xp})}: diagonal, predictor-specific, but identical across observations, \deqn{Qx_i = \mathrm{diag}\left(1/\sigma_{x1}^2, \dots, 1/\sigma_{xp}^2\right) \text{ for all } i}
\eqn{p \times p} matrix \eqn{\Sigma_x}: global covariance allowing correlated predictor errors, identical across observations, \deqn{Qx_i = \Sigma_x^{-1} \text{ for all } i}
\eqn{n \times p} matrix (or, when \eqn{p = 1}, a length-\eqn{n} vector as shorthand): observation-specific standard deviations \eqn{\sigma_{x,i1}, \dots, \sigma_{x,ip}}, giving a genuinely observation-varying, diagonal-only precision \deqn{Qx_i = \mathrm{diag}\left(1/\sigma_{x,i1}^2, \dots, 1/\sigma_{x,ip}^2\right)}
Cross-predictor correlation is not supported at the per-observation level (only globally, via the \eqn{p \times p} form above), since that would require a separate \eqn{p \times p} matrix per observation.
In every case except the last, \eqn{Qx_i} is the same fixed matrix/vector for all observations and is simply written \eqn{Qx}. 
If \eqn{\sigma_x} is (nearly) zero for a predictor that is in fact measured without error (for example a time index), supplying a small \code{sigma_x} makes the fit approach ordinary nonlinear least squares in that direction. With the default unit precisions, the relative scale of the predictors and the response directly determines how much of the misfit is attributed to \eqn{x} rather than \eqn{y}.
  
\bold{Response precision (\code{sigma_y}, \code{weights}):} The baseline response precision is \eqn{1/\sigma_{y,i}^2}, using either the common \code{sigma_y} or, if a length-\eqn{n} vector is supplied, the observation-specific value. If observation weights \eqn{w_i} are supplied, \deqn{Qyy_i = w_i / \sigma_{y,i}^2} so that observation weights act multiplicatively on the response precision and thus enter the objective \eqn{S} directly (affecting both the foot-point corrections and the parameter estimates), rather than acting solely as conventional \code{\link{nls}}-style regression weights. Weights act on the response side only; the predictor precision \eqn{Qx_i} is not scaled by them.

\bold{Relation to ODRPACK weights (\code{WE}, \code{WD}).} ODRPACK takes a response weight \code{WE} and a predictor weight \code{WD} for every observation, both as precisions (inverse variances). In \code{onls} they correspond to \eqn{Qyy_i = WE_i} and \eqn{Qx_i = WD_i} (a diagonal matrix if \eqn{p > 1}), and are specified as \code{weights = WE} (together with the default \code{sigma_y = 1}; alternatively \code{sigma_y = 1/sqrt(WE)} without \code{weights}) and \code{sigma_x = 1/sqrt(WD)} (for \eqn{p > 1} an \eqn{n \times p}{n x p} matrix with entries \eqn{1/\sqrt{WD_{ij}}}{1/sqrt(WD_ij)}). ODRPACK always scales the standard errors by the residual variance, which corresponds to \code{known_sigma = FALSE} (supplying \code{sigma_x} or \code{sigma_y} otherwise sets \code{known_sigma = TRUE}). A full weight matrix that differs between observations (a three-dimensional \code{WD} in ODRPACK) is not supported, only a global covariance matrix (see \code{sigma_x}). The demo \code{gauge} (\code{\link{doGauge}}) fits the ODRPACK example with separate \code{WE} and \code{WD} weights.

\bold{Algorithm.} 
    
1) \emph{Warm start.} A conventional nonlinear least-squares model is fitted using \code{\link[minpack.lm]{nlsLM}} (optionally weighted by \code{weights}; \code{fixed} parameters are substituted into the model formula) to obtain starting values \eqn{\theta^{(0)}}. If this fit fails, a warning is issued and the raw \code{start} values are used instead.\cr
2) \emph{Joint Levenberg-Marquardt problem.} With \eqn{L_i} the factor satisfying \eqn{L_i^T L_i = Qx_i} (elementwise square root for a diagonal \eqn{Qx_i}, Cholesky factor for a correlated global \eqn{\Sigma_x}), define the residual vector of length \eqn{n + np} 
\deqn{r(\theta, \delta) = \left( \left\{ Qyy_i^{1/2}\left(y_i - f(x_i + \delta_i, \theta)\right) \right\}_{i=1}^{n}, \; \left\{ L_i \delta_i \right\}_{i=1}^{n} \right)}
so that \eqn{S = r^T r}. Starting at \eqn{(\theta^{(0)}, \delta = 0)}, the Levenberg-Marquardt iteration (\enc{Moré}{More} 1978; \code{\link[minpack.lm]{nls.lm}}) repeatedly solves \deqn{\left(J^T J + \mu D^2\right) \Delta z = -J^T r \qquad z = (\theta_{free}, \delta_1, \dots, \delta_n)} with a trust-region controlled damping \eqn{\mu} and a diagonal scaling \eqn{D} (the running maximum of the Jacobian column norms, or \code{control$diag}). The Jacobian of \eqn{r} has the sparse 'arrow' structure of ODRPACK,
\deqn{J = \left( \begin{array}{cc} -W_y^{1/2} F_\theta & -W_y^{1/2} G \\ 0 & L_x \end{array} \right)}{J = [ -Wy^(1/2) F_theta, -Wy^(1/2) G ; 0, Lx ]} 
where \eqn{W_y = \mathrm{diag}(Qyy_i)}, \eqn{F_\theta} is the \eqn{n \times q_{free}} matrix with rows \eqn{\partial f(\xi_i,\theta)/\partial\theta^T}, \eqn{G} is the \eqn{n \times np} block-diagonal matrix whose \eqn{i}-th row holds \eqn{g_i^T = \partial f(\xi_i,\theta)/\partial x^T} in the columns belonging to \eqn{\delta_i}, and \eqn{L_x = \mathrm{blockdiag}(L_1,\dots,L_n)}. Each correction \eqn{\delta_i} thus affects only observation \eqn{i}. The two algorithms differ only in how the step is computed:\cr
\emph{\code{algorithm = "structured"}} (default for \eqn{n > 100} or with parameter bounds). In the normal equations, the \eqn{\delta}-part is block diagonal, with one \eqn{p \times p}{p x p} block \eqn{D_i = Qyy_i g_i g_i^T + Qx_i} (+ damping) per observation. These blocks are eliminated (Schur complement, as in ODRPACK), so that each step needs one \eqn{q_{free} \times q_{free}}{q x q} solve and \eqn{n} small block solves (Sherman-Morrison for diagonal \eqn{Qx_i}, a Cholesky factorization vectorized over the observations for a correlated \eqn{\Sigma_x}):
\deqn{\Big(A + \mu D_\theta^2 - \sum_i B_i \tilde D_i^{-1} B_i^T\Big)\Delta\theta = -g_\theta + \sum_i B_i \tilde D_i^{-1} g_{\delta,i} \qquad \Delta\delta_i = -\tilde D_i^{-1}\big(g_{\delta,i} + B_i^T \Delta\theta\big)}{(A + mu D^2 - sum B_i Dt_i^-1 B_i') dtheta = -g_theta + sum B_i Dt_i^-1 g_delta_i,  ddelta_i = -Dt_i^-1 (g_delta_i + B_i' dtheta),}
with \eqn{A = F_\theta^T W_y F_\theta}, \eqn{B_i = Qyy_i \, \partial f_i/\partial\theta \, g_i^T}, \eqn{\tilde D_i = D_i + \mu D_{\delta,i}^2}, and \eqn{g_\theta}, \eqn{g_{\delta,i}} the parts of \eqn{J^T r}. The damping \eqn{\mu} is chosen by MINPACK's trust-region rule (\enc{Moré}{More} 1978): \eqn{\mu = 0} (a Gauss-Newton step) when that step lies inside the trust region, otherwise \eqn{\mu} such that \eqn{\|D \Delta z\|} matches its radius. Parameter bounds \code{lower}/\code{upper} are handled by holding a parameter whose step would leave its bounds at the bound and re-solving for the others. The cost per iteration is \eqn{O(n p^3 + n q^2 + q^3)}, linear in \eqn{n}; a fit with \eqn{n = 100\,000} takes seconds. For small \eqn{n} its per-iteration overhead in \R is higher than that of the compiled \pkg{minpack.lm} code, so near-degenerate fits that run thousands of iterations (for example some bootstrap resamples in \code{\link{confint.onls}}) can take longer than with \code{"dense"}.\cr
\emph{\code{algorithm = "dense"}}. The joint problem is passed as one nonlinear least-squares problem with the full \eqn{(n + np) \times (q + np)}{(n + np) x (q + np)} Jacobian to \code{\link[minpack.lm]{nls.lm}}, whose dense QR factorization costs \eqn{O((n + q)^3)} per iteration (about 15 s per fit at \eqn{n = 1000}). It is kept for comparison and is the only algorithm that supports the foot-point bounds \code{extend}/\code{window}.\cr
Both algorithms give the same estimates up to the convergence tolerance on well-posed problems. For models with several local minima (for example periodic models with free periods) they can end in different minima, because their iterates differ slightly.\cr
3) \emph{Analytic derivatives.} The blocks \eqn{F_\theta} and \eqn{G} are supplied to the solver exactly, by symbolic differentiation of the model formula (\code{\link{deriv}}) when possible, and by central finite differences otherwise (for example if the formula calls a user-defined function). This matters: MINPACK's own forward-difference Jacobian uses a step proportional to the size of each unknown, which becomes smaller than the floating-point resolution of \eqn{f(x_i + \delta_i)} once a correction \eqn{\delta_i} is small, and would then drive that correction to exactly zero, i.e. to a foot point that is not orthogonal.\cr
4) \emph{Iterations and convergence.} The structured algorithm iterates until convergence (codes 1 to 3, as in \code{nls.lm}, giving \code{convergence = 0}) or until \code{control$outer_max} iterations are used (\code{convergence = 1}, with a warning). \pkg{minpack.lm} caps a single \code{nls.lm} call at 1024 iterations, so the dense algorithm is restarted from its last iterate (which also resets the trust region) until it converges (\code{nls.lm} codes 1 to 4), a restart no longer reduces the objective, or the budget is exhausted. Convergence is declared with the tolerances \code{ftol} and \code{ptol} (default \code{1e-10}). If the joint optimization fails to run at all, the NLS estimates with zero corrections are returned, with a warning.

\bold{Stationarity (first-order) conditions.} At a solution, \eqn{\partial S/\partial \delta_i = 0} and \eqn{\partial S/\partial \theta = 0} give
\deqn{Qx_i (\hat\xi_i - x_i) = Qyy_i \left(y_i - f(\hat\xi_i, \hat\theta)\right) \nabla_x f(\hat\xi_i, \hat\theta) \qquad \sum_{i=1}^{n} Qyy_i \left(y_i - f(\hat\xi_i, \hat\theta)\right) \frac{\partial f(\hat\xi_i, \hat\theta)}{\partial \theta} = 0}
For unit precisions the first condition states that the vector from the foot point to the observation is orthogonal to the model surface at the foot point, which is what \code{\link{check_o}} verifies, as the angle between the residual vector and the model surface after each axis has been rescaled by its precision.

\bold{Optional foot-point bounds (\code{extend}, \code{window}).} By default the corrections \eqn{\delta_i} are unbounded, as in ODRPACK. For single-predictor models, supplying \code{extend} confines the foot points to \eqn{[\text{XLOW}, \text{XUPP}]},
\deqn{\text{XLOW} = \min(x) - \text{extend[1]} \times \mathrm{range}(x) \qquad \text{XUPP} = \max(x) + \text{extend[2]} \times \mathrm{range}(x)}
by bounding \eqn{\delta_i}. If, in addition, \code{window} is given and \eqn{n > 25}, the foot point of observation \eqn{i} (in sorted predictor order) is confined to \eqn{[x_{(i-\text{window}+1)},\, x_{(i+\text{window}-1)}]}, clipped at the two ends of the data to \code{XLOW} and \code{XUPP}. With bounds, the bounded variant of the Levenberg-Marquardt algorithm is used, and a foot point sitting at a bound is generally not orthogonal. The bounded variant can also stop before the optimum, in particular with \code{window}; \code{onls} then warns when foot points are neither orthogonal nor at their bound. \code{window} has no effect unless \code{extend} is supplied; for multivariate models (\eqn{p>1}) neither argument is used. These bounds are rarely needed and are kept for backward compatibility and for models whose foot points must be kept on one branch of the curve.

\bold{Fixed parameters:} If \code{fixed} is supplied, parameters flagged \code{TRUE} are held, throughout both the NLS warm start and the joint optimization, at their user-supplied \code{start} value (not at the NLS warm-start estimate). Fixed parameters contribute neither to the residual degrees of freedom (\code{df_resid} \eqn{= n - q_{free}}, with \eqn{q_{free}} the number of free parameters) nor to the Jacobian columns used for \code{vcov}/\code{std_errors}; their standard error is reported as \code{0}. Observations with zero weight do not count towards \code{df_resid}.

\bold{Covariance of the parameter estimates.} The covariance is the ODRPACK covariance: the \eqn{\theta}-block of the inverse Gauss-Newton matrix \eqn{(J^T J)^{-1}} of the joint problem (second-derivative terms are neglected, as in ODRPACK). Eliminating the \eqn{\delta}-blocks of \eqn{J^T J} exactly (Schur complement, using the Woodbury identity) yields the compact form used by \code{onls}, valid for any number of predictors \eqn{p}: 
\deqn{\widehat{\mathrm{Var}}(\hat\theta_{free}) = \left(F_\theta^T W F_\theta\right)^{-1} \qquad \widehat{\mathrm{Var}}(\hat\theta_{free}) = \hat\sigma^2 \left(F_\theta^T W F_\theta\right)^{-1} \quad}
for known_sigma = \code{TRUE} and known_sigma = \code{FALSE}, respectively, and with \eqn{F_\theta} evaluated at the final foot points \eqn{\hat\xi_i} and estimates \eqn{\hat\theta}, \eqn{W = \mathrm{diag}(w_1,\dots,w_n)} and the effective (\dQuote{effective-variance}) observation weights
\deqn{w_i = \left(Qyy_i^{-1} + g_i^T \, Qx_i^{-1} \, g_i\right)^{-1} \qquad g_i = \nabla_x f(\hat\xi_i, \hat\theta)}
which combine the response error and the predictor error propagated through the local model gradient (for \eqn{p = 1}, \eqn{w_i = (1/Qyy_i + s_i^2/Qx_i)^{-1}} with slope \eqn{s_i}). 
Here \eqn{\hat\sigma^2 = \left(\sum_i d_i^2\right)/\text{df\_resid}} is the reduced chi-square (\code{reduced_chisq} in 'Value'). Standard errors are the square roots of the diagonal.

\bold{Starting values, local minima and degenerate solutions.} The objective can have several local minima (for example for periodic models with free periods, or for models that contain a limiting case such as the Richards curve, whose limit is the Gompertz curve), and the solution found depends on \code{start}. Starting values that make the model insensitive to some parameters (for example a sigmoid whose exponent underflows for all observations, so that the model is a constant) create a degenerate stationary point at which the solver may report convergence; \code{onls} then issues a warning naming the parameters that have no measurable influence on the fit, and the results should be discarded. A low orthogonal residual sum of squares is not by itself evidence of a good fit: with unit precisions and a steep model, an ODR fit can lower the objective by shifting observations horizontally. Compare the vertical residuals (\code{\link{residuals}}) with those of the NLS fit at \code{parNLS}, and choose \code{sigma_x}/\code{sigma_y} to reflect the actual measurement errors.

The returned object contains both classical vertical residual information and orthogonal-distance quantities such as foot points (\code{xi}, \code{yi}; the predictor corrections are \code{xi - X_obs}) and weighted orthogonal distances.

\bold{IMPORTANT:} If not all points are orthogonal to the fitted curve, \code{\link{print.onls}} gives a \dQuote{FAILED: Only X out of Y fitted points are orthogonal} message. In this case, it is suggested to conduct a more detailed analysis using \code{\link{check_o}}. Because orthogonality of the foot points is the stationarity condition \eqn{\partial S/\partial \delta_i = 0}, the most common cause is a convergence tolerance that is too loose for observations with very small residuals, for which the orthogonality angle is very sensitive to the foot-point position: tightening \code{control = list(ftol = 1e-12, ptol = 1e-12)} (or increasing \code{control$outer_max} if the iteration budget was exhausted) will normally resolve it. Foot points held at a bound set through \code{extend}/\code{window} are not expected to be orthogonal.

\bold{ALSO IMPORTANT:} Regular \R-like notation such as \code{lm(y ~ x1 + x2)} or \code{lm(y ~ x1 * x2)} needs to be changed to classical notation \eqn{y = b_0 + b_1x1 + b_2x2} or \eqn{y = b_0 + b_1x1 + b_2x2 + b_3x1x2}, respectively.

The resulting orthogonal model houses information in respect to the (classical) vertical residuals as well as the (minimized, precision-weighted) orthogonal residuals.

The following functions use the \bold{vertical} residuals:\cr
\code{\link{deviance}}\cr
\code{\link{fitted}}\cr
\code{\link{residuals}}\cr
\code{\link{logLik}}

The following functions use the \bold{orthogonal} residuals:\cr
\code{\link{deviance_o}}\cr
\code{\link{residuals_o}}\cr
\code{\link{logLik_o}}
\code{\link{anova.onls}} 
}

\value{
  An orthogonal fit of class \code{onls} (also inheriting from \code{nls}), a list with the following elements:\cr
  \item{\code{call}, \code{data}, \code{formula}}{The matched call, the \code{data} argument as given (unevaluated) and the model formula.}
  \item{\code{algorithm}}{The algorithm used, \code{"structured"} or \code{"dense"} (see \code{algorithm}).}
  \item{\code{coefficients}}{Named vector of the final parameter estimates.}
  \item{\code{std_errors}, \code{vcov}}{Standard errors (\code{0} for \code{fixed} parameters) and covariance matrix of the estimates; see 'Details'.}
  \item{\code{parNLS}}{Parameter estimates of the initial nonlinear least-squares fit (the starting values, if the warm start failed; fixed parameters at their starting values).}
  \item{\code{convergence}, \code{convInfo}}{Convergence code (\code{0} if the joint optimization converged, termination codes 1 to 4; \code{1} otherwise) and details: \code{isConv}, \code{finIter} (total number of Levenberg-Marquardt iterations), \code{finTol} (the \code{ftol} used), \code{stopCode} and \code{stopMessage}.}
  \item{\code{objective}}{The minimized objective \eqn{S = \sum_i d_i^2}.}
  \item{\code{reduced_chisq}}{Reduced chi-square \eqn{S/\mathrm{df\_resid}}.}
  \item{\code{known_sigma}}{Logical: were the variances treated as known?}
  \item{\code{df_resid}}{Residual degrees of freedom: the number of observations with non-zero weight minus the number of free parameters.}
  \item{\code{n}}{Number of observations used (after \code{subset} and \code{na.action}).}
  \item{\code{X_obs}, \code{y_obs}}{Observed predictors (\eqn{n \times p}{n x p} matrix) and responses.}
  \item{\code{xi}, \code{yi}}{Foot points \eqn{\hat\xi_i} (\eqn{n \times p}{n x p} matrix) and foot-point responses \eqn{f(\hat\xi_i, \hat\theta)}; see \code{\link{x0}}.}
  \item{\code{fittedONLS}}{Model evaluated at the \emph{observed} predictors, \eqn{f(x_i, \hat\theta)}; the vertical residuals are \code{y_obs - fittedONLS} (\code{\link{residuals}}).}
  \item{\code{orth_dist}}{Weighted orthogonal distances \eqn{d_i}.}
  \item{\code{weights}}{Observation weights.}
  \item{\code{Q_yy}, \code{Lx}}{Response precisions \eqn{Qyy_i} (length \eqn{n}) and predictor precisions as an \eqn{n \times p \times p}{n x p x p} array of upper Cholesky factors, \eqn{L_i^T L_i = Qx_i}: one representation for every form of \code{sigma_x}.}
  \item{\code{sigma_x}, \code{sigma_y}, \code{per_obs_x}}{\code{sigma_x} and \code{sigma_y} as used (after \code{subset}/\code{na.action}), and whether \code{sigma_x} was given per observation.}
  \item{\code{fixed}, \code{lower}, \code{upper}, \code{control}, \code{extend}, \code{window}}{The arguments of the same names, used when the model is refitted by \code{\link{profile.onls}} and \code{\link{confint.onls}}.}
  \item{\code{na.action}, \code{dataClasses}}{Information on omitted observations, and the classes of the model variables (used by \code{\link{predict}}).}
  \item{\code{ortho}}{Logical vector, one value per observation: the \code{Ortho} column of \code{\link{check_o}} (\code{NA} if it could not be computed).}
All per-observation elements are in the original observation order. Version 0.3-0 removed the duplicated elements of earlier versions (\code{parONLS}, \code{x0}/\code{y0}, \code{resid_o}, \code{pred}/\code{resp}, \code{residONLS}, \code{residNLS}, \code{fittedNLS}, \code{grad}, \code{QR}, \code{Q_x}, \code{Sigma_x}, \code{model}, ...); use \code{coef}, \code{\link{x0}}/\code{\link{y0}}, \code{\link{residuals_o}}, \code{\link{residuals}} and \code{\link{fitted}} instead.
}

\author{
  Andrej-Nikolai Spiess
}   

\references{
A stable and efficient algorithm for nonlinear orthogonal distance regression.\cr
Boggs PT, Byrd RH and Schnabel RB.\cr
SIAM J Sci Stat Comput (1987), \bold{8}: 1052-1078.\cr
\doi{10.1137/0908085}.

The Levenberg-Marquardt algorithm: implementation and theory.\cr
More JJ.\cr
In: Watson GA (ed.), Numerical Analysis, Lecture Notes in Mathematics \bold{630}, Springer (1978): 105-116.\cr
\doi{10.1007/BFb0067700}.

Orthogonal Distance Regression.\cr
Boggs PT and Rogers JE.\cr
NISTIR (1990), \bold{89-4197}: 1-15.\cr
\url{https://static.scipy.org/doc/external/odr_ams.pdf}.
  
User's Reference Guide for ODRPACK Version 2.01\cr
Software for Weighted Orthogonal Distance Regression.\cr
Boggs PT, Byrd RH, Rogers JE and Schnabel RB.\cr
NISTIR (1992), \bold{4834}: 1-113.\cr
\url{https://static.scipy.org/doc/external/odrpack_guide.pdf}.

ALGORITHM 676 ODRPACK: Software for Weighted Orthogonal Distance Regression.\cr
Boggs PT, Donaldson JR, Byrd RH and Schnabel RB.\cr
ACM Trans Math Soft (1989), \bold{15}, 348-364.\cr
\doi{10.1145/76909.76913}.

Omnivariant generalized least squares regression.\cr
Daeron M and Vermeesch P, Chemical Geology (2024), \bold{647}: 121881.\cr
\doi{10.1016/j.chemgeo.2023.121881}.

An analysis of the Total Least Squares problem.\cr
Golub GH and Van Loan CF.\cr
SIAM Journal on Numerical Analysis (1980), \bold{17}, 883-893.\cr
\doi{10.1137/0717073}.

Least squares fitting of a straight line with correlated errors.\cr
York D.\cr
Earth and Planetary Science Letters (1968), \bold{5}, 320-324.\cr
\doi{10.1016/S0012-821X(68)80059-7}.
}

\examples{
## 1. The DNase data from 'nls', use all generic functions.
DNase1 <- subset(DNase, Run == 1)
set.seed(123)
DNase1$density <- sapply(DNase1$density, function(x) rnorm(1, x, 0.1 * x))
mod1 <- onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
             data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1))
print(mod1)
plot(mod1)
summary(mod1)
predict(mod1, newdata = data.frame(conc = 6))
logLik(mod1)
deviance(mod1)
formula(mod1)
weights(mod1)
df.residual(mod1)
fitted(mod1)
residuals(mod1)
vcov(mod1)
coef(mod1)

## 2a. Update model
DNase2 <- DNase1
DNase2$conc <- DNase2$conc * 2
mod2a <- update(mod1, data = DNase2)
print(mod2a)

## 2b. Example with a fixed parameter
## => Asym = 3.
mod2b <- onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
             data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1), 
             fixed = c(TRUE, FALSE, FALSE))
print(mod2b)

## 3. Multivariate example: matched curvature, 
## low noise, decorrelated predictors 
set.seed(123)
n  <- 25
x1 <- runif(n, 1, 5)     
x2 <- runif(n, 1, 5)
b1_true <- 5
b2_true <- 2
b3_true <- 1.5
z_true <- b1_true + b2_true * x1 + b3_true * x2^2
z  <- z_true + rnorm(n, 0, 2)      
x1 <- x1 + rnorm(n, 0, 0.5)
x2 <- x2 + rnorm(n, 0, 0.5)
DAT <- data.frame(x1 = x1, x2 = x2, z = z)
mod3 <- onls(z ~ b1 + b2 * x1 + b3 * x2^2, data = DAT, 
             start = list(b1 = 1, b2 = 1, b3 = 1), trace = TRUE)
print(mod3)

## The reference tests against published ODR results, closed-form solutions, ODRPACK and
## the NIST certified values are collected in one demo: doGauge(), or
## demo("gauge", package = "onls").

\donttest{
## 4. The two algorithms: same fit, very different cost for large n.
set.seed(1)
n <- 400
xt <- seq(0, 10, length.out = n)
DAT <- data.frame(x = xt + rnorm(n, 0, 0.2), y = 5/(1 + exp(-(xt - 5))) + rnorm(n, 0, 0.2))
system.time(m_s <- onls(y ~ b1/(1 + exp(-(x - b2)/b3)), data = DAT,
                        start = list(b1 = 4, b2 = 4, b3 = 1.5)))   # auto: structured for n > 100
system.time(m_d <- onls(y ~ b1/(1 + exp(-(x - b2)/b3)), data = DAT,
                        start = list(b1 = 4, b2 = 4, b3 = 1.5), algorithm = "dense"))
rbind(structured = coef(m_s), dense = coef(m_d))
}
}   

\keyword{optimize}
\keyword{models}
\keyword{nonlinear}
