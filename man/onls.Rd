\name{onls}
\alias{onls}
\encoding{latin1}

\title{Weighted Orthogonal Nonlinear Least-Squares Regression (Orthogonal Distance Regression)}

\description{
  Fits nonlinear orthogonal-distance regression (ODR, also called errors-in-variables) models in the manner of ODRPACK (Boggs, Byrd, Rogers and Schnabel; 1987): the model parameters and one foot-point correction per observation are estimated \emph{simultaneously} by minimizing the same weighted sum of squared response and predictor residuals that ODRPACK minimizes, using a Levenberg-Marquardt algorithm on the joint problem. 
  In contrast to ordinary nonlinear least squares, \code{onls()} allows measurement error in both predictors and response variables. 
  Predictor uncertainty may be specified through predictor-specific standard deviations, observation-specific error structures, or full covariance matrices allowing correlated predictor errors. Response uncertainty may be supplied globally or on a per-observation basis and can be combined with observation weights. 
  The algorithm supports single- and multi-predictor nonlinear models, parameter bounds, fixed parameters, and covariance estimation for the fitted parameters that is equivalent to the ODRPACK (Gauss-Newton) covariance.
}

\usage{
  onls(formula, data, start = NULL, weights = NULL, sigma_x = NULL, sigma_y = 1, 
       known_sigma = NULL, extend = NULL, window = NULL, control = list(), 
       lower = NULL, upper = NULL, fixed = NULL, subset = NULL, na.action = NULL, 
       trace = FALSE)
}

\arguments{
  \item{formula}{A two-sided nonlinear model \code{\link{formula}} describing the response as a nonlinear function of one or more predictors and unknown model parameters. The response must be a plain variable name (no transformation of the left-hand side).}
  \item{data}{A data frame containing the variables appearing in \code{formula}.}
  \item{start}{A named list (or named vector) of starting values for all model parameters.}
  \item{weights}{Optional non-negative observation weights. The weights are incorporated into the orthogonal-distance criterion through the response precision matrix and therefore affect both the foot-point corrections and the parameter estimates. When supplied, \code{weights} must have length equal to the number of observations.}
  \item{sigma_x}{Specification of predictor measurement error, as a standard deviation (not a variance or a weight). Accepts \code{NULL}, a scalar, a length-\eqn{p} vector, a \eqn{p \times p} covariance matrix, or an \eqn{n \times p} matrix of observation-specific standard deviations (a length-\eqn{n} vector is accepted as shorthand for the latter when \eqn{p = 1}). See 'Details'.}
  \item{sigma_y}{Response measurement error, as a standard deviation. May be either a single positive value applied to all observations or a length-\eqn{n} vector of observation-specific standard deviations.}
  \item{known_sigma}{Logical indicating whether the supplied (or default) \code{sigma_x}/\code{sigma_y} should be treated as known measurement-error standard deviations rather than as estimated from the data. If \code{NULL} (default), it is set to \code{TRUE} if \code{sigma_x} is supplied or \code{sigma_y} is explicitly passed by the user (regardless of its value), and \code{FALSE} otherwise. Standard errors are scaled by the reduced chi-square statistic (see \code{reduced_chisq} in 'Value') only when \code{known_sigma = FALSE}; when \code{known_sigma = TRUE} the supplied variances are assumed to already be on the correct absolute scale and no such scaling is applied.}
  \item{extend}{Optional numeric vector of length one or two. If \code{NULL} (default), the foot-point corrections are unbounded, as in ODRPACK. If supplied (single-predictor models only), the foot points are confined to a range extended beyond the observed predictor range by the given fractions (lower, upper; a single value is used for both ends). See 'Details'.}
  \item{window}{Optional integer window width, used only if \code{extend} is supplied, for single-predictor models with \eqn{n > 25}: restricts the foot point of each observation to the neighborhood of the corresponding sorted observations. See 'Details'.}
  \item{control}{Optional list of optimization settings. Recognized elements are \code{outer_max}, the \emph{total} budget of Levenberg-Marquardt iterations (default \code{5000}, spent in restarts of at most 1024 iterations each, see 'Details'); \code{ftol} and \code{ptol}, the relative-reduction and relative-step convergence tolerances passed to \code{\link[minpack.lm]{nls.lm.control}} (default \code{1e-10} each); \code{diag}, an optional length-\eqn{q} vector of positive parameter scaling factors (foot-point corrections are given the scale 1; if omitted, MINPACK's adaptive scaling is used); and \code{lower}/\code{upper}, alternatives to the arguments of the same name.}
  \item{lower}{Optional vector of lower bounds for the model parameters, of length \eqn{q} (same ordering as \code{start}, fixed parameters included).}
  \item{upper}{Optional vector of upper bounds for the model parameters, of length \eqn{q} (same ordering as \code{start}, fixed parameters included).}
  \item{fixed}{Optional logical vector with the same length and ordering as \code{start}. Parameters marked \code{TRUE} are held fixed at their supplied starting values throughout optimization. At least one parameter must be free.}
  \item{subset}{Optional specification of a subset of observations to be used for fitting. Row-aligned \code{weights}, \code{sigma_x} and \code{sigma_y} are subset accordingly.}
  \item{na.action}{Function indicating how missing values should be handled. Row-aligned \code{weights}, \code{sigma_x} and \code{sigma_y} follow the removal of observations.}
  \item{trace}{Logical. If \code{TRUE}, progress messages from the NLS warm start and the orthogonal-distance optimization (including restarts) are printed.}
}

\details{
\bold{Model and objective.} Assume a nonlinear model \eqn{y = f(x, \theta)}, where \eqn{x} is a vector of \eqn{p} predictors and \eqn{\theta} is a vector of \eqn{q} unknown model parameters.
Unlike ordinary nonlinear least squares, which assumes error only in the response, ODR assumes that both the response and the predictors are observed with error. 
For observation \eqn{i}, let \eqn{x_i} be the observed predictor vector, \eqn{y_i} the observed response, \eqn{\delta_i} the unknown correction to the predictors (so that \eqn{\xi_i = x_i + \delta_i} is the \emph{foot point} on the model surface),
\eqn{Qyy_i} the response precision and \eqn{Qx_i} the predictor precision matrix (both constructed below). 
The parameters and all \eqn{n} corrections are estimated jointly by minimizing
\deqn{S(\theta, \delta_1, \dots, \delta_n) = \sum_{i=1}^{n} \left[ Qyy_i \left(y_i - f(x_i + \delta_i, \theta)\right)^2 + \delta_i^T Qx_i \, \delta_i \right].}
This is the explicit ODR problem of ODRPACK, with \eqn{Qyy_i} playing the role of the response weights (\code{WE}) and \eqn{Qx_i} the role of the predictor weights (\code{WD}). It has \eqn{q_{free} + np} unknowns, where \eqn{q_{free}} is the number of free (non-\code{fixed}) parameters.
For fixed \eqn{\theta}, minimizing \eqn{S} over \eqn{\delta_i} alone gives the weighted squared orthogonal distance of observation \eqn{i} to the model surface,
\deqn{d_i^2 = Qyy_i \left[y_i - f(\hat\xi_i, \theta) \right]^2 + (\hat\xi_i - x_i)^T Qx_i (\hat\xi_i - x_i),}
so that the estimate \eqn{\hat\theta} is the minimizer of \eqn{\sum_i d_i^2}; the \eqn{d_i} are returned as \code{orth_dist} and \code{resid_o}, and \eqn{\sum_i d_i^2} as \code{objective}.
  
\bold{Predictor precision (\code{sigma_x}):}\cr
The predictor precision \eqn{Qx_i} used above is constructed from \code{sigma_x}, which is always supplied on the standard-deviation scale:\cr
\code{NULL}: unit variance for every predictor, \eqn{Qx_i = I_p} for all \eqn{i}\cr
scalar \eqn{\sigma}: isotropic, \eqn{Qx_i = I_p/\sigma^2} for all \eqn{i}\cr
length-\eqn{p} vector \eqn{(\sigma_{x1}, \dots, \sigma_{xp})}: diagonal, predictor-specific, but identical across observations, \deqn{Qx_i = \mathrm{diag}\left(1/\sigma_{x1}^2, \dots, 1/\sigma_{xp}^2\right) \text{ for all } i.}
\eqn{p \times p} matrix \eqn{\Sigma_x}: global covariance allowing correlated predictor errors, identical across observations, \deqn{Qx_i = \Sigma_x^{-1} \text{ for all } i.}
\eqn{n \times p} matrix (or, when \eqn{p = 1}, a length-\eqn{n} vector as shorthand): observation-specific standard deviations \eqn{\sigma_{x,i1}, \dots, \sigma_{x,ip}}, giving a genuinely observation-varying, diagonal-only precision \deqn{Qx_i = \mathrm{diag}\left(1/\sigma_{x,i1}^2, \dots, 1/\sigma_{x,ip}^2\right).}
Cross-predictor correlation is not supported at the per-observation level (only globally, via the \eqn{p \times p} form above), since that would require a separate \eqn{p \times p} matrix per observation.
In every case except the last, \eqn{Qx_i} is the same fixed matrix/vector for all observations and is simply written \eqn{Qx}. 
If \eqn{\sigma_x} is (nearly) zero for a predictor that is in fact measured without error (for example a time index), supplying a small \code{sigma_x} makes the fit approach ordinary nonlinear least squares in that direction. With the default unit precisions, the relative scale of the predictors and the response directly determines how much of the misfit is attributed to \eqn{x} rather than \eqn{y}.
  
\bold{Response precision (\code{sigma_y}, \code{weights}):} The baseline response precision is \eqn{1/\sigma_{y,i}^2}, using either the common \code{sigma_y} or, if a length-\eqn{n} vector is supplied, the observation-specific value. If observation weights \eqn{w_i} are supplied, \deqn{Qyy_i = w_i / \sigma_{y,i}^2,} so that observation weights act multiplicatively on the response precision and thus enter the objective \eqn{S} directly (affecting both the foot-point corrections and the parameter estimates), rather than acting solely as conventional \code{\link{nls}}-style regression weights. Weights act on the response side only; the predictor precision \eqn{Qx_i} is not scaled by them.

\bold{Relation to ODRPACK weights (\code{WE}, \code{WD}).} ODRPACK takes a response weight \code{WE} and a predictor weight \code{WD} for every observation, both as precisions (inverse variances). In \code{onls} they correspond to \eqn{Qyy_i = WE_i} and \eqn{Qx_i = WD_i} (a diagonal matrix if \eqn{p > 1}), and are specified as \code{weights = WE} (together with the default \code{sigma_y = 1}; alternatively \code{sigma_y = 1/sqrt(WE)} without \code{weights}) and \code{sigma_x = 1/sqrt(WD)} (for \eqn{p > 1} an \eqn{n \times p}{n x p} matrix with entries \eqn{1/\sqrt{WD_{ij}}}{1/sqrt(WD_ij)}). ODRPACK always scales the standard errors by the residual variance, which corresponds to \code{known_sigma = FALSE} (supplying \code{sigma_x} or \code{sigma_y} otherwise sets \code{known_sigma = TRUE}). A full weight matrix that differs between observations (a three-dimensional \code{WD} in ODRPACK) is not supported, only a global covariance matrix (see \code{sigma_x}). See example 12.

\bold{Algorithm.} 
    
1) \emph{Warm start.} A conventional nonlinear least-squares model is fitted using \code{\link[minpack.lm]{nlsLM}} (optionally weighted by \code{weights}; \code{fixed} parameters are substituted into the model formula) to obtain starting values \eqn{\theta^{(0)}}. If this fit fails, a warning is issued and the raw \code{start} values are used instead.\cr
2) \emph{Joint Levenberg-Marquardt problem.} With \eqn{L_i} the factor satisfying \eqn{L_i^T L_i = Qx_i} (elementwise square root for a diagonal \eqn{Qx_i}, Cholesky factor for a correlated global \eqn{\Sigma_x}), define the residual vector of length \eqn{n + np} 
\deqn{r(\theta, \delta) = \left( \left\{ Qyy_i^{1/2}\left(y_i - f(x_i + \delta_i, \theta)\right) \right\}_{i=1}^{n}, \; \left\{ L_i \delta_i \right\}_{i=1}^{n} \right),}
so that \eqn{S = r^T r}. Starting at \eqn{(\theta^{(0)}, \delta = 0)}, the Levenberg-Marquardt iteration (\enc{Moré}{More} 1978; \code{\link[minpack.lm]{nls.lm}}) repeatedly solves \deqn{\left(J^T J + \mu D^2\right) \Delta z = -J^T r, \qquad z = (\theta_{free}, \delta_1, \dots, \delta_n),} with a trust-region controlled damping \eqn{\mu} and scaling \eqn{D}. The Jacobian of \eqn{r} has the sparse 'arrow' structure of ODRPACK,
\deqn{J = \left( \begin{array}{cc} -W_y^{1/2} F_\theta & -W_y^{1/2} G \\ 0 & L_x \end{array} \right),}{J = [ -Wy^(1/2) F_theta, -Wy^(1/2) G ; 0, Lx ],} 
where \eqn{W_y = \mathrm{diag}(Qyy_i)}, \eqn{F_\theta} is the \eqn{n \times q_{free}} matrix with rows \eqn{\partial f(\xi_i,\theta)/\partial\theta^T}, \eqn{G} is the \eqn{n \times np} block-diagonal matrix whose \eqn{i}-th row holds \eqn{g_i^T = \partial f(\xi_i,\theta)/\partial x^T} in the columns belonging to \eqn{\delta_i}, and \eqn{L_x = \mathrm{blockdiag}(L_1,\dots,L_n)}. Each correction \eqn{\delta_i} thus affects only observation \eqn{i}. ODRPACK exploits this structure to make each step cheap; \code{onls} delegates the step to the general dense solver of \pkg{minpack.lm}, so the cost per iteration grows faster with \eqn{n} than in ODRPACK.\cr
3) \emph{Analytic derivatives.} The blocks \eqn{F_\theta} and \eqn{G} are supplied to the solver exactly, by symbolic differentiation of the model formula (\code{\link{deriv}}) when possible, and by central finite differences otherwise (for example if the formula calls a user-defined function). This matters: MINPACK's own forward-difference Jacobian uses a step proportional to the size of each unknown, which becomes smaller than the floating-point resolution of \eqn{f(x_i + \delta_i)} once a correction \eqn{\delta_i} is small, and would then drive that correction to exactly zero, i.e. to a foot point that is not orthogonal.\cr
4) \emph{Restarts and convergence.} \pkg{minpack.lm} caps a single \code{nls.lm} call at 1024 iterations. The solver is restarted from its last iterate (which also resets the trust region) until it converges (\code{nls.lm} codes 1 to 4, giving \code{convergence = 0}), a restart no longer reduces the objective, or the budget is exhausted (\code{convergence = 1}, with a warning). Convergence is declared with the tolerances \code{ftol} and \code{ptol} (default \code{1e-10}). If the joint optimization fails to run at all, the NLS estimates with zero corrections are returned, with a warning.

\bold{Stationarity (first-order) conditions.} At a solution, \eqn{\partial S/\partial \delta_i = 0} and \eqn{\partial S/\partial \theta = 0} give
\deqn{Qx_i (\hat\xi_i - x_i) = Qyy_i \left(y_i - f(\hat\xi_i, \hat\theta)\right) \nabla_x f(\hat\xi_i, \hat\theta), \qquad \sum_{i=1}^{n} Qyy_i \left(y_i - f(\hat\xi_i, \hat\theta)\right) \frac{\partial f(\hat\xi_i, \hat\theta)}{\partial \theta} = 0.}
For unit precisions the first condition states that the vector from the foot point to the observation is orthogonal to the model surface at the foot point, which is what \code{\link{check_o}} verifies (as an angle or, for weighted fits, through the relative residual of this equation).

\bold{Optional foot-point bounds (\code{extend}, \code{window}).} By default the corrections \eqn{\delta_i} are unbounded, as in ODRPACK. For single-predictor models, supplying \code{extend} confines the foot points to \eqn{[\text{XLOW}, \text{XUPP}]},
\deqn{\text{XLOW} = \min(x) - \text{extend[1]} \times \mathrm{range}(x), \qquad \text{XUPP} = \max(x) + \text{extend[2]} \times \mathrm{range}(x),}
by bounding \eqn{\delta_i}. If, in addition, \code{window} is given and \eqn{n > 25}, the foot point of observation \eqn{i} (in sorted predictor order) is confined to \eqn{[x_{(i-\text{window}+1)},\, x_{(i+\text{window}-1)}]}, clipped at the two ends of the data to \code{XLOW} and \code{XUPP}. With bounds, the bounded variant of the Levenberg-Marquardt algorithm is used, and a foot point sitting at a bound is generally not orthogonal. \code{window} has no effect unless \code{extend} is supplied; for multivariate models (\eqn{p>1}) neither argument is used. These bounds are rarely needed and are kept for backward compatibility and for models whose foot points must be kept on one branch of the curve.

\bold{Fixed parameters:} If \code{fixed} is supplied, parameters flagged \code{TRUE} are held, throughout both the NLS warm start and the joint optimization, at their user-supplied \code{start} value (not at the NLS warm-start estimate). Fixed parameters contribute neither to the residual degrees of freedom (\code{df_resid} \eqn{= n - q_{free}}, with \eqn{q_{free}} the number of free parameters) nor to the Jacobian columns used for \code{vcov}/\code{std_errors}; their standard error is reported as \code{0}.

\bold{Covariance of the parameter estimates.} The covariance is the ODRPACK covariance: the \eqn{\theta}-block of the inverse Gauss-Newton matrix \eqn{(J^T J)^{-1}} of the joint problem (second-derivative terms are neglected, as in ODRPACK). Eliminating the \eqn{\delta}-blocks of \eqn{J^T J} exactly (Schur complement, using the Woodbury identity) yields the compact form used by \code{onls}, valid for any number of predictors \eqn{p}: 
\deqn{\widehat{\mathrm{Var}}(\hat\theta_{free}) = \left(F_\theta^T W F_\theta\right)^{-1} \qquad \widehat{\mathrm{Var}}(\hat\theta_{free}) = \hat\sigma^2 \left(F_\theta^T W F_\theta\right)^{-1} \quad}
for known_sigma = \code{TRUE} and known_sigma = \code{FALSE}, respectively, and with \eqn{F_\theta} evaluated at the final foot points \eqn{\hat\xi_i} and estimates \eqn{\hat\theta}, \eqn{W = \mathrm{diag}(w_1,\dots,w_n)} and the effective (\dQuote{effective-variance}) observation weights
\deqn{w_i = \left(Qyy_i^{-1} + g_i^T \, Qx_i^{-1} \, g_i\right)^{-1}, \qquad g_i = \nabla_x f(\hat\xi_i, \hat\theta),}
which combine the response error and the predictor error propagated through the local model gradient (for \eqn{p = 1}, \eqn{w_i = (1/Qyy_i + s_i^2/Qx_i)^{-1}} with slope \eqn{s_i}). 
Here \eqn{\hat\sigma^2 = \left(\sum_i d_i^2\right)/\text{df\_resid}} is the reduced chi-square (\code{reduced_chisq} in 'Value'). Standard errors are the square roots of the diagonal.

\bold{Starting values, local minima and degenerate solutions.} The objective can have several local minima (for example for periodic models with free periods, or for models that contain a limiting case such as the Richards curve, whose limit is the Gompertz curve), and the solution found depends on \code{start}. Starting values that make the model insensitive to some parameters (for example a sigmoid whose exponent underflows for all observations, so that the model is a constant) create a degenerate stationary point at which the solver may report convergence; \code{onls} then issues a warning naming the parameters that have no measurable influence on the fit, and the results should be discarded. A low orthogonal residual sum of squares is not by itself evidence of a good fit: with unit precisions and a steep model, an ODR fit can lower the objective by shifting observations horizontally. Compare the vertical residuals (\code{residONLS}) with those of the NLS fit (\code{residNLS}), and choose \code{sigma_x}/\code{sigma_y} to reflect the actual measurement errors.

The returned object contains both classical vertical residual information and orthogonal-distance quantities such as foot points, predictor corrections and weighted orthogonal distances.

\bold{IMPORTANT:} If not all points are orthogonal to the fitted curve, \code{\link{print.onls}} gives a \dQuote{FAILED: Only X out of Y fitted points are orthogonal} message. In this case, it is suggested to conduct a more detailed analysis using \code{\link{check_o}}. Because orthogonality of the foot points is the stationarity condition \eqn{\partial S/\partial \delta_i = 0}, the most common cause is a convergence tolerance that is too loose for observations with very small residuals, for which the orthogonality angle is very sensitive to the foot-point position: tightening \code{control = list(ftol = 1e-12, ptol = 1e-12)} (or increasing \code{control$outer_max} if the iteration budget was exhausted) will normally resolve it. Foot points held at a bound set through \code{extend}/\code{window} are not expected to be orthogonal.

\bold{ALSO IMPORTANT:} Regular \R-like notation such as \code{lm(y ~ x1 + x2)} or \code{lm(y ~ x1 * x2)} needs to be changed to classical notation \eqn{y = b_0 + b_1x1 + b_2x2} or \eqn{y = b_0 + b_1x1 + b_2x2 + b_3x1x2}, respectively.

The resulting orthogonal model houses information in respect to the (classical) vertical residuals as well as the (minimized Euclidean) orthogonal residuals.

The following functions use the \bold{vertical} residuals:\cr
\code{\link{deviance}}\cr
\code{\link{fitted}}\cr
\code{\link{residuals}}\cr
\code{\link{logLik}}

The following functions use the \bold{orthogonal} residuals:\cr
\code{\link{deviance_o}}\cr
\code{\link{residuals_o}}\cr
\code{\link{logLik_o}}
}

\value{
  An orthogonal fit of class \code{onls} with the following list items:\cr
  \item{\code{data}}{Original \code{data} argument.}
  \item{\code{call}}{Matched function call.}
  \item{\code{convInfo}}{Convergence information from the joint optimization: \code{isConv} (logical), \code{finIter} (total number of Levenberg-Marquardt iterations over all restarts), \code{finTol} (the \code{ftol} used), \code{stopCode} and \code{stopMessage} (termination code and message of the last \code{\link[minpack.lm]{nls.lm}} call).}
  \item{\code{na.action}}{Information on omitted observations.}
  \item{\code{dataClasses}}{Classes of model variables.}
  \item{\code{model}}{Model frame used for fitting.}
  \item{\code{formula}}{Model formula.}
  \item{\code{parNLS}}{Parameter estimates from the initial nonlinear least-squares fit (the starting values, if the warm start failed; fixed parameters at their starting values).}
  \item{\code{parONLS}}{Final orthogonal-distance parameter estimates.}
  \item{\code{x0}, \code{y0}}{Coordinates \eqn{\hat\xi_i, f(\hat\xi_i,\hat\theta)} of the fitted foot points on the model curve/surface, in sorted predictor order (see 'Details'; for \eqn{p>1}, sorting is a no-op and the original order is used).}
  \item{\code{fittedONLS}}{Model evaluated at the \emph{observed} predictors, \eqn{f(x_i, \hat\theta)}, using the final ODR parameter estimates.}
  \item{\code{fittedNLS}}{Model evaluated at the observed predictors using the initial NLS warm-start estimates.}
  \item{\code{residONLS}}{Vertical residuals \eqn{y_i - f(x_i,\hat\theta)} from the orthogonal model (observed \eqn{x}, not the foot point).}
  \item{\code{residNLS}}{Vertical residuals from the initial nonlinear least-squares fit.}
  \item{\code{resid_o}}{Weighted orthogonal distances \eqn{d_i} (see 'Details') at the final estimates, in sorted predictor order. Numerically equivalent to \code{orth_dist}, but reordered; compare to \code{orth_dist} for the original observation order.}
  \item{\code{pred}, \code{resp}}{Predictor and response values used internally during fitting (sorted order for \eqn{p=1}).}
  \item{\code{grad}}{Jacobian \eqn{\partial f(\hat\xi_i,\theta)/\partial\theta} of the model with respect to \emph{all} model parameters (fixed and free), evaluated at the final foot points \eqn{\hat\xi_i} and estimates \eqn{\hat\theta} (not at the observed predictors), in the original observation order. Only its free-parameter columns are used to construct \code{vcov}.}
  \item{\code{QR}}{QR decomposition of \code{grad}.}
  \item{\code{weights}}{Observation weights used in fitting, in original observation order.}
  \item{\code{control}}{The resolved control settings actually passed to \code{\link[minpack.lm]{nls.lm}} (via \code{\link[minpack.lm]{nls.lm.control}}) for each restart, not the raw \code{control} argument.}
  \item{\code{coefficients}}{Named vector of final parameter estimates.}
  \item{\code{std_errors}}{Approximate standard errors of parameter estimates (\code{0} for \code{fixed} parameters); see 'Details' for the covariance used.}
  \item{\code{vcov}}{Estimated parameter covariance matrix; see 'Details'.}
  \item{\code{xi}}{Estimated foot-point coordinates \eqn{\hat\xi_i} in the original observation order.}
  \item{\code{x_corrections}}{Estimated predictor corrections \eqn{\hat\delta_i = \hat\xi_i - x_i}.}
  \item{\code{orth_dist}}{Final weighted orthogonal distances \eqn{d_i}, in original observation order (see \code{resid_o} for the sorted-order equivalent).}
  \item{\code{objective}}{Minimized objective value \eqn{S(\hat\theta, \hat\delta) = \sum_i d_i^2}.}
  \item{\code{Q_x}}{Predictor precision structure used in fitting.}
  \item{\code{Sigma_x}}{Predictor covariance structure corresponding to \code{Q_x}.}
  \item{\code{Q_yy}}{Response precision vector.}
  \item{\code{is_diag_x}}{Logical indicating whether predictor precision is diagonal.}
  \item{\code{per_obs_x}}{Logical indicating whether predictor precision varies by observation.}
  \item{\code{sigma_x}}{User-supplied predictor error specification (after \code{subset}/\code{na.action}).}
  \item{\code{sigma_y}}{User-supplied response error specification (after \code{subset}/\code{na.action}).}
  \item{\code{known_sigma}}{Logical indicating whether measurement variances were treated as known.}
  \item{\code{reduced_chisq}}{Reduced chi-square statistic based on the orthogonal-distance criterion.}
  \item{\code{pred_names}}{Predictor variable names.}
  \item{\code{param_names}}{Model parameter names.}
  \item{\code{start}}{Starting parameter values.}
  \item{\code{fixed}}{Logical vector identifying fixed parameters.}
  \item{\code{lower}, \code{upper}}{Parameter bounds supplied to the optimizer.}
  \item{\code{y_obs}, \code{X_obs}}{Observed response and predictor values.}
  \item{\code{resp_name}}{Response variable name.}
  \item{\code{n}, \code{p}, \code{q}}{Numbers of observations, predictors and model parameters.}
  \item{\code{df_resid}}{Residual degrees of freedom, \eqn{n} minus the number of \emph{free} (non-\code{fixed}) parameters.}
  \item{\code{convergence}}{Convergence code, \code{0} if the joint optimization converged (\code{nls.lm} termination codes 1 to 4), \code{1} otherwise.}
  \item{\code{iters}}{Total number of Levenberg-Marquardt iterations, summed over all restarts (\code{NA} if the joint optimization failed to run).}
  \item{\code{type}}{Model type, currently \code{"nonlinear"}.}
  \item{\code{ortho}}{Orthogonality diagnostics returned by \code{\link{check_o}} (\code{NA} if they could not be computed).}
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

\donttest{
## Reference tests comparing to pivotal literature
## 4. Example from odrpack_guide.pdf, 2.C.i, pages 39ff.
x <- c(0, 0, 5, 7, 7.5, 10, 16, 26, 30, 34, 34.5, 100)
y <- c(1265, 1263.6, 1258, 1254, 1253, 1249.8, 1237, 1218, 1220.6, 
1213.8, 1215.5, 1212)
DAT <- data.frame(x, y)
mod4 <- onls(y ~ b1 + b2 * (exp(b3 * x) -1)^2, data = DAT, 
             start = list(b1 = 1500, b2 = -50, b3 = -0.1))
deviance_o(mod4) # 21.445 as on page 47
summary(mod4) # 1264.65481 (1.03492) / -54.01838 (1.583992) / -0.08785 (6.33222E-3) as on page 48

## 5. Example from Algorithm 676: ODRPACK, page 355 + 356.
x <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 85, 90, 95, 100, 105)
y <- c(4.14, 8.52, 16.31, 32.18, 64.62, 98.76, 151.13, 224.74, 341.35, 
       423.36, 522.78, 674.32, 782.04, 920.01)
DAT <- data.frame(x, y)
mod5 <- onls(y ~ b1 * 10^(b2 * x/(b3 + x)), data = DAT, 
             start = list(b1 = 1, b2 = 5, b3 = 100))
deviance_o(mod5) # 15.263 as on page 363
summary(mod5) # 4.4879 (0.56876) / 7.1882 (0.69504) / 221.8383 (37.2313) as on page 363

## 6. Example with bounds from simple_example.f90
## in https://www.netlib.org/toms/869.zip.
x <- c(0.982, 1.998, 4.978, 6.01)
y <- c(2.7, 7.4, 148.0, 403.0)
DAT <- data.frame(x, y)
mod6 <- onls(y ~ b1 * exp(b2 * x), data = DAT, 
            start = list(b1 = 2, b2 = 0.5), 
            lower = c(0, 0), upper = c(10, 0.9))
coef(mod6) # 1.4376 / 0.9 ## Different to reference 1.6334 / 0.9
deviance_o(mod6) # 0.1919 => lower RSS than original ODRPACK with 0.2674!

## 7. Example similar to Deming regression
## Comparison to XLstat
## https://help.xlstat.com/6650-run-deming-regression-compare-methods-excel
x <- c(9.8, 9.7, 10.7, 10.9, 12.4, 12.5, 12.8, 12.8, 12.9, 13.3, 
       13.4, 13.5, 13.7, 14.9, 15.2, 15.5)
y <- c(10.1, 11.4, 10.8, 11.3, 11.8, 12.1, 12.3, 13.6, 14.2, 14.4,
       14.6, 15.3, 15.5, 15.8, 16.2, 16.5)
DAT <- data.frame(x, y)
mod7 <- onls(y ~ a + b * x, data = DAT, start = list(a = 2, b = 3))
print(mod7) ## -1.909 / 1.208 as on webpage
plot(mod7)

## 8. Linear multivariate model, using the closed-form Total Least Squares 
## (TLS) solution from Golub & Van Loan (1980)
tls_fit <- function(X, y) {
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  yc <- y - mean(y)
  xbar <- colMeans(X); ybar <- mean(y)
  Z <- cbind(Xc, yc)                
  SVD <- svd(Z)
  v <- SVD$v[, p + 1L]               
  v_x <- v[1:p]; v_y <- v[p + 1L]
  slope <- -v_x / v_y; intercept <- ybar - sum(slope * xbar)
    list(intercept = intercept, slope = setNames(slope, colnames(X)),
       singular_values = SVD$d)
}

set.seed(11)
n <- 40
x1_true <- runif(n, 0, 10)
x2_true <- runif(n, 0, 10)
b0_true <- 3; b1_true <- 1.5; b2_true <- -0.8
y_true  <- b0_true + b1_true * x1_true + b2_true * x2_true
x1 <- x1_true + rnorm(n, 0, 0.5)
x2 <- x2_true + rnorm(n, 0, 0.5)
y  <- y_true  + rnorm(n, 0, 0.5)
DAT <- data.frame(x1 = x1, x2 = x2, y = y)
TLS <- tls_fit(DAT[, c("x1", "x2")], DAT$y)
mod8 <- onls(y ~ b0 + b1 * x1 + b2 * x2, data = DAT, start = list(b0 = 1, b1 = 1, b2 = 1))
TLS_vec  <- c(b0 = TLS$intercept, b1 = TLS$slope[["x1"]], b2 = TLS$slope[["x2"]])
ONLS_vec <- coef(mod8)[c("b0", "b1", "b2")]
print(data.frame(TLS_closed_form = TLS_vec, onls = ONLS_vec, 
      abs_diff = abs(TLS_vec - ONLS_vec))) # all equal

## 9. Pearson (1901) / York (1966) -> "Pearson's data with York's weights"
#   Intercept: 5.47991  (SE 0.29497)
#   Slope:    -0.48053  (SE 0.05799)
x <- c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
sd_x <- 1/sqrt(c(1000.0, 1000.0, 500.0, 800.0, 200.0, 80.0, 60.0, 20.0, 1.8, 1.0))
sd_y <- 1/sqrt(c(1.0, 1.8, 4.0, 8.0, 20.0, 20.0, 70.0, 70.0, 100.0, 500.0))
DAT <- data.frame(x = x, y = y)
mod9 <- onls(y ~ b0 + b1*x, data = DAT,
                    start = list(b0 = 5, b1 = -0.5),
                    sigma_x = sd_x, sigma_y = sd_y)
summary(mod9) # 5.47991 (0.29497) / -0.48053 (0.05799) as in paper

## 10. Daeron & Vermeesch (2024), Table 3 / Figure 2C toy example.
x <- c(9, 19, 31, 41)
y <- c(21, 31, 39, 49)
DAT <- data.frame(x = x, y = y)
mod10 <- onls(y ~ a + b * x, data = DAT, start = list(a = 10, b = 1), sigma_x = 1, sigma_y = 1)
summary(mod10) # 13.71 / 0.851 as in Table 3 of paper 

## 11. Full predictor covariance (correlated predictor errors), compared to the closed-form 
## generalized Total Least Squares (TLS) solution. Whitening the predictors with the Cholesky 
## factor L of the precision matrix (t(L) \%*\% L = solve(Sigma)) and scaling y by 1/sigma_y 
## turns the problem into plain TLS, which is solved by an SVD.
set.seed(123)
n <- 40
Sigma <- matrix(c(0.25, 0.15, 0.15, 0.16), 2)  # correlation of predictor errors = 0.75
sigma_y <- 0.3
xt <- cbind(runif(n, 0, 10), runif(n, 0, 10))
E <- matrix(rnorm(2 * n), n) \%*\% chol(Sigma)
DAT <- data.frame(x1 = xt[, 1] + E[, 1], x2 = xt[, 2] + E[, 2],
                  y = 3 + 1.5 * xt[, 1] - 0.8 * xt[, 2] + rnorm(n, 0, sigma_y))
mod12 <- onls(y ~ b0 + b1 * x1 + b2 * x2, data = DAT, 
              start = list(b0 = 1, b1 = 1, b2 = 1), sigma_x = Sigma, sigma_y = sigma_y,
              control = list(ftol = 1e-13, ptol = 1e-13))
L <- chol(solve(Sigma))
U <- as.matrix(DAT[, c("x1", "x2")]) \%*\% t(L) # whitened predictors
v <- DAT$y / sigma_y
Z <- cbind(scale(U, scale = FALSE), v - mean(v))
V <- svd(Z)$v[, 3]
w <- -V[1:2] / V[3]
gTLS <- c(sigma_y * (mean(v) - sum(w * colMeans(U))), sigma_y * drop(t(L) \%*\% w))
print(data.frame(gen_TLS = gTLS, onls = coef(mod12), 
                 abs_diff = abs(gTLS - coef(mod12)))) # all equal
  
## 12. ODRPACK's separate weights WE (response) and WD (predictor).
## ODRPACK takes both as precisions (inverse variances), per observation. In onls(), WE is 
## passed as 'weights' (with the default sigma_y = 1) and WD through sigma_x = 1/sqrt(WD) 
## (for p > 1: an n x p matrix 1/sqrt(WD)). known_sigma = FALSE gives ODRPACK's scaling of 
## the standard errors by the residual variance.
set.seed(123)
n <- 30
xt <- seq(0.5, 10, length.out = n)
WD <- runif(n, 0.5, 4)   # predictor weights
WE <- runif(n, 0.5, 4)   # response weights
x <- xt + rnorm(n, 0, 0.3/sqrt(WD))
y <- 2 * exp(-0.3 * xt) + 0.5 + rnorm(n, 0, 0.03/sqrt(WE))
DAT <- data.frame(x, y)
mod12 <- onls(y ~ a * exp(-b * x) + c, data = DAT, start = list(a = 1.5, b = 0.2, c = 0.3),
              weights = WE, sigma_x = 1/sqrt(WD), known_sigma = FALSE)
summary(mod12) 
}
}   

\keyword{optimize}
\keyword{models}
\keyword{nonlinear}
