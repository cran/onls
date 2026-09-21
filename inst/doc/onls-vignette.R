## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.width = 6,
  fig.height = 5.5
)

## ----setup--------------------------------------------------------------------
library(onls)

## ----univariate-fit-----------------------------------------------------------
DNase1 <- subset(DNase, Run == 1)
set.seed(1)
DNase1$density <- sapply(DNase1$density, function(x) rnorm(1, x, 0.1 * x))

mod_uni <- onls(density ~ Asym / (1 + exp((xmid - log(conc)) / scal)),
                data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1))
print(mod_uni)

## ----univariate-summary-------------------------------------------------------
summary(mod_uni)

## ----univariate-plot, fig.alt = "Orthogonal vs vertical fit with foot-point segments"----
plot(mod_uni)

## ----univariate-check---------------------------------------------------------
check_o(mod_uni, plot = FALSE)

## ----foot-points--------------------------------------------------------------
FP <- data.frame(x = mod_uni$pred, x0 = x0(mod_uni), y = mod_uni$resp, y0 = y0(mod_uni))
FP$dist <- sqrt((FP$x - FP$x0)^2 + (FP$y - FP$y0)^2)
head(FP)

## for unit precisions, the sum of squared distances is the minimized objective
all.equal(sum(FP$dist^2), deviance_o(mod_uni))

## ----orthogonal-residuals-----------------------------------------------------
d_o <- residuals_o(mod_uni)
all.equal(sort(d_o), sort(FP$dist), check.attributes = FALSE)
all.equal(sum(d_o^2), deviance_o(mod_uni))

## side by side with the vertical residuals of the same fit
head(data.frame(vertical = residuals(mod_uni), orthogonal = d_o))

## ----orthogonal-residuals-na--------------------------------------------------
DNase1_na <- DNase1
DNase1_na$density[c(3, 10)] <- NA
mod_na <- onls(density ~ Asym / (1 + exp((xmid - log(conc)) / scal)),
               data = DNase1_na, start = list(Asym = 3, xmid = 0, scal = 1),
               na.action = na.exclude)
residuals_o(mod_na)

## ----plot-quadratic, fig.alt = "Noisy quadratic with orthogonal and vertical fit and foot-point segments"----
set.seed(123)
x <- 1:20
y <- 10 + 3 * x^2 + rnorm(20, 0, 50)
DAT_quad <- data.frame(x, y)
mod_quad <- onls(y ~ a + b * x^2, data = DAT_quad, start = list(a = 10, b = 3))
plot(mod_quad)

## ----plot-zoom, fig.alt = "Zoom into the left half of the quadratic fit, without the NLS curve"----
plot(mod_quad, fitted.nls = FALSE, xlim = c(0, 10), asp = FALSE)

## ----weighted-fit-------------------------------------------------------------
x <- c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
sd_x <- 1 / sqrt(c(1000, 1000, 500, 800, 200, 80, 60, 20, 1.8, 1.0))
sd_y <- 1 / sqrt(c(1.0, 1.8, 4.0, 8.0, 20, 20, 70, 70, 100, 500))
DAT_py <- data.frame(x = x, y = y)

mod_w <- onls(y ~ b0 + b1 * x, data = DAT_py, start = list(b0 = 5, b1 = -0.5),
              sigma_x = sd_x, sigma_y = sd_y)
summary(mod_w)   # intercept 5.480 (0.295), slope -0.481 (0.058), matching York's published values

## ----weighted-check-----------------------------------------------------------
check_o(mod_w, plot = FALSE)

## ----we-wd--------------------------------------------------------------------
set.seed(7)
n <- 30
xt <- seq(0.5, 10, length.out = n)
WD <- runif(n, 0.5, 4)   # predictor weights
WE <- runif(n, 0.5, 4)   # response weights
x <- xt + rnorm(n, 0, 0.3 / sqrt(WD))
y <- 2 * exp(-0.3 * xt) + 0.5 + rnorm(n, 0, 0.03 / sqrt(WE))
DAT_we <- data.frame(x, y)

mod_we <- onls(y ~ a * exp(-b * x) + c, data = DAT_we,
               start = list(a = 1.5, b = 0.2, c = 0.3),
               weights = WE, sigma_x = 1 / sqrt(WD), known_sigma = FALSE)
summary(mod_we)   # 1.98131 (0.03347) / 0.28984 (0.01468) / 0.48480 (0.02986), as scipy.odr with we = WE, wd = WD

## ----we-wd-objective----------------------------------------------------------
f <- function(x, b) b[1] * exp(-b[2] * x) + b[3]
xi <- mod_we$xi[, 1]
all.equal(sum(WE * (y - f(xi, coef(mod_we)))^2 + WD * (xi - x)^2), mod_we$objective)

## ----odrpack-guide------------------------------------------------------------
x <- c(0, 0, 5, 7, 7.5, 10, 16, 26, 30, 34, 34.5, 100)
y <- c(1265, 1263.6, 1258, 1254, 1253, 1249.8, 1237, 1218, 1220.6, 1213.8, 1215.5, 1212)
DAT_guide <- data.frame(x, y)

mod_guide <- onls(y ~ b1 + b2 * (exp(b3 * x) - 1)^2, data = DAT_guide,
                  start = list(b1 = 1500, b2 = -50, b3 = -0.1))
deviance_o(mod_guide)   # 21.445, as on page 47 of the guide
summary(mod_guide)      # 1264.65481 (1.03492) / -54.01838 (1.583992) / -0.08785 (6.33222E-3), as on page 48

## ----odrpack-guide-check------------------------------------------------------
check_o(mod_guide, plot = FALSE)

## ----algorithm-676------------------------------------------------------------
x <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 85, 90, 95, 100, 105)
y <- c(4.14, 8.52, 16.31, 32.18, 64.62, 98.76, 151.13, 224.74, 341.35,
       423.36, 522.78, 674.32, 782.04, 920.01)
DAT_676 <- data.frame(x, y)

mod_676 <- onls(y ~ b1 * 10^(b2 * x / (b3 + x)), data = DAT_676,
                start = list(b1 = 1, b2 = 5, b3 = 100))
deviance_o(mod_676)   # 15.263, as on page 363
summary(mod_676)      # 4.4879 (0.56876) / 7.1882 (0.69504) / 221.8383 (37.2313), as on page 363

## ----daeron-------------------------------------------------------------------
DAT_dv <- data.frame(x = c(9, 19, 31, 41), y = c(21, 31, 39, 49))
mod_dv <- onls(y ~ a + b * x, data = DAT_dv, start = list(a = 10, b = 1),
               sigma_x = 1, sigma_y = 1)
summary(mod_dv)   # 13.71 / 0.8516 (the exact TLS slope); Table 3 of the paper lists 13.71 / 0.851

## ----bivariate-fit------------------------------------------------------------
set.seed(2024)
n  <- 60
x1 <- runif(n, -5, 5)
x2 <- runif(n, -5, 5)

a_true <- 3; b_true <- 2; c_true <- 4
z <- c_true * sqrt(1 + (x1 / a_true)^2 + (x2 / b_true)^2) + rnorm(n, 0, 0.3)

x1 <- x1 + rnorm(n, 0, 0.2)
x2 <- x2 + rnorm(n, 0, 0.15)
DAT_hyp <- data.frame(x1 = x1, x2 = x2, z = z)

mod_hyp <- onls(z ~ c * sqrt(1 + (x1 / a)^2 + (x2 / b)^2), data = DAT_hyp,
                start = list(a = 2, b = 2, c = 3),
                sigma_x = c(0.2, 0.15), sigma_y = 0.3)
summary(mod_hyp)   # expect a, b, c close to 3, 2, 4

## ----bivariate-check----------------------------------------------------------
check_o(mod_hyp, plot = FALSE)

## ----bivariate-plot, eval = FALSE---------------------------------------------
# plot(mod_hyp)

## ----half-dome----------------------------------------------------------------
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
DAT_dome <- data.frame(x1 = x1, x2 = x2, z = z)
maxrad <- max(sqrt(x1^2 + x2^2))

mod_dome <- onls(z ~ sqrt(r^2 - x1^2 - x2^2), data = DAT_dome,
                 start = list(r = r_true),
                 sigma_x = c(0.1, 0.1), sigma_y = 0.15,
                 lower = maxrad * 1.05, upper = 100)
summary(mod_dome)   # r close to 6
check_o(mod_dome, plot = FALSE)   # all orthogonal to the dome surface

## ----half-dome-plot, eval = FALSE---------------------------------------------
# plot(mod_dome)   # renders as a visibly round dome

## ----classical-lm-style, eval = FALSE-----------------------------------------
# ## This is NOT valid onls() syntax -- shown only for comparison.
# lm(z ~ x1 + I(x2^2) + sqrt(x3) + log(x4 + 1))

## ----classical-onls-style, eval = FALSE---------------------------------------
# z ~ b0 + b1 * x1 + b2 * x2^2 + b3 * sqrt(x3) + b4 * log(x4 + 1)

## ----multivariate-fit---------------------------------------------------------
set.seed(99)
n  <- 60
x1 <- runif(n, 0, 10)
x2 <- runif(n, 0, 5)
x3 <- runif(n, 2, 10)   # kept away from 0: sqrt() needs non-negative arguments
x4 <- runif(n, 2, 10)   # kept away from -1: log(x4 + 1) needs x4 + 1 > 0

b0 <- 2; b1 <- 0.8; b2 <- 0.5; b3 <- 2; b4 <- 3
z <- b0 + b1 * x1 + b2 * x2^2 + b3 * sqrt(x3) + b4 * log(x4 + 1) + rnorm(n, 0, 0.5)

sd_x <- c(0.3, 0.2, 0.3, 0.3)
x1 <- x1 + rnorm(n, 0, sd_x[1])
x2 <- x2 + rnorm(n, 0, sd_x[2])
x3 <- x3 + rnorm(n, 0, sd_x[3])
x4 <- x4 + rnorm(n, 0, sd_x[4])
DAT_mv <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, z = z)

mod_mv <- onls(z ~ b0 + b1 * x1 + b2 * x2^2 + b3 * sqrt(x3) + b4 * log(x4 + 1),
               data = DAT_mv,
               start = list(b0 = 1, b1 = 1, b2 = 1, b3 = 1, b4 = 1),
               sigma_x = sd_x, sigma_y = 0.5)
summary(mod_mv)   # expect b0..b4 close to 2, 0.8, 0.5, 2, 3

## ----multivariate-check-------------------------------------------------------
check_o(mod_mv, plot = FALSE)

## ----multivariate-grid, fig.width = 7, fig.height = 7, fig.alt = "Grid of four partial-dependence panels"----
plot(mod_mv)

## ----multivariate-panel, fig.alt = "Single full-size panel for one predictor"----
plot(mod_mv, panel = "x2")

## ----full-covariance----------------------------------------------------------
set.seed(2026)
n <- 40
Sigma <- matrix(c(0.25, 0.15, 0.15, 0.16), 2)     # correlation of predictor errors = 0.75
sigma_y <- 0.3
xt <- cbind(runif(n, 0, 10), runif(n, 0, 10))
E <- matrix(rnorm(2 * n), n) %*% chol(Sigma)
DAT_cov <- data.frame(x1 = xt[, 1] + E[, 1], x2 = xt[, 2] + E[, 2],
                      y = 3 + 1.5 * xt[, 1] - 0.8 * xt[, 2] + rnorm(n, 0, sigma_y))

mod_cov <- onls(y ~ b0 + b1 * x1 + b2 * x2, data = DAT_cov,
                start = list(b0 = 1, b1 = 1, b2 = 1), sigma_x = Sigma, sigma_y = sigma_y,
                control = list(ftol = 1e-13, ptol = 1e-13))

## closed-form generalized TLS
L <- chol(solve(Sigma))
U <- as.matrix(DAT_cov[, c("x1", "x2")]) %*% t(L)      # whitened predictors
v <- DAT_cov$y / sigma_y
Z <- cbind(scale(U, scale = FALSE), v - mean(v))
V <- svd(Z)$v[, 3]
w <- -V[1:2] / V[3]
gTLS <- c(sigma_y * (mean(v) - sum(w * colMeans(U))), sigma_y * drop(t(L) %*% w))
print(data.frame(gen_TLS = gTLS, onls = coef(mod_cov),
                 abs_diff = abs(gTLS - coef(mod_cov))))   # all equal

## ----full-covariance-diag-----------------------------------------------------
mod_cov_d <- onls(y ~ b0 + b1 * x1 + b2 * x2, data = DAT_cov,
                  start = list(b0 = 1, b1 = 1, b2 = 1), sigma_x = sqrt(diag(Sigma)),
                  sigma_y = sigma_y)
coef(mod_cov_d)

## ----deming-------------------------------------------------------------------
x <- c(9.8, 9.7, 10.7, 10.9, 12.4, 12.5, 12.8, 12.8, 12.9, 13.3,
       13.4, 13.5, 13.7, 14.9, 15.2, 15.5)
y <- c(10.1, 11.4, 10.8, 11.3, 11.8, 12.1, 12.3, 13.6, 14.2, 14.4,
       14.6, 15.3, 15.5, 15.8, 16.2, 16.5)
DAT_dem <- data.frame(x, y)

mod_dem <- onls(y ~ a + b * x, data = DAT_dem, start = list(a = 2, b = 3))
print(mod_dem)   # -1.909 / 1.208 as on the webpage

## ----tls----------------------------------------------------------------------
tls_fit <- function(X, y) {
  X <- as.matrix(X)
  p <- ncol(X)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  yc <- y - mean(y)
  xbar <- colMeans(X); ybar <- mean(y)
  SVD <- svd(cbind(Xc, yc))
  v <- SVD$v[, p + 1L]
  slope <- -v[1:p] / v[p + 1L]
  list(intercept = ybar - sum(slope * xbar), slope = setNames(slope, colnames(X)))
}

set.seed(11)
n <- 40
x1_true <- runif(n, 0, 10)
x2_true <- runif(n, 0, 10)
y_true  <- 3 + 1.5 * x1_true - 0.8 * x2_true
DAT_tls <- data.frame(x1 = x1_true + rnorm(n, 0, 0.5),
                      x2 = x2_true + rnorm(n, 0, 0.5),
                      y  = y_true  + rnorm(n, 0, 0.5))

TLS <- tls_fit(DAT_tls[, c("x1", "x2")], DAT_tls$y)
mod_tls <- onls(y ~ b0 + b1 * x1 + b2 * x2, data = DAT_tls,
                start = list(b0 = 1, b1 = 1, b2 = 1))
TLS_vec  <- c(b0 = TLS$intercept, b1 = TLS$slope[["x1"]], b2 = TLS$slope[["x2"]])
ONLS_vec <- coef(mod_tls)[c("b0", "b1", "b2")]
print(data.frame(TLS_closed_form = TLS_vec, onls = ONLS_vec,
                 abs_diff = abs(TLS_vec - ONLS_vec)))   # equal to solver tolerance

## ----bounds-------------------------------------------------------------------
DAT_bnd <- data.frame(x = c(0.982, 1.998, 4.978, 6.01),
                      y = c(2.7, 7.4, 148.0, 403.0))
mod_bnd <- onls(y ~ b1 * exp(b2 * x), data = DAT_bnd,
                start = list(b1 = 2, b2 = 0.5),
                lower = c(0, 0), upper = c(10, 0.9))
coef(mod_bnd)          # 1.4376 / 0.9, different to the reference 1.6334 / 0.9
deviance_o(mod_bnd)    # 0.1919, lower than the 0.2674 of the original ODRPACK

## ----fixed--------------------------------------------------------------------
mod_fix <- onls(density ~ Asym / (1 + exp((xmid - log(conc)) / scal)),
                data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1),
                fixed = c(TRUE, FALSE, FALSE))
print(mod_fix)

## ----control------------------------------------------------------------------
mod_ctrl <- onls(y ~ b1 + b2 * (exp(b3 * x) - 1)^2, data = DAT_guide,
                 start = list(b1 = 1500, b2 = -50, b3 = -0.1),
                 control = list(ftol = 1e-12, ptol = 1e-12, outer_max = 2000))
coef(mod_ctrl)
mod_ctrl$convInfo$isConv
mod_ctrl$convInfo$finIter   # total Levenberg-Marquardt iterations

## ----degenerate---------------------------------------------------------------
x <- 1:15
y <- c(16.08, 33.83, 65.80, 97.20, 191.55, 326.20, 386.87, 520.53,
       590.03, 651.92, 724.93, 699.56, 689.96, 637.56, 717.41)
DAT_rich <- data.frame(x, y)

mod_flat <- withCallingHandlers(
  onls(y ~ b1 / (1 + exp(b2 - b3 * x))^(1 / b4), data = DAT_rich,
       start = list(b1 = 10, b2 = -1, b3 = 7, b4 = 9)),
  warning = function(w) {
    message("Warning: ", conditionMessage(w))
    invokeRestart("muffleWarning")
  })

## ----gompertz-----------------------------------------------------------------
mod_gomp <- onls(y ~ b1 * exp(-exp(c - b3 * x)), data = DAT_rich,
                 start = list(b1 = 750, c = 2, b3 = 0.5))
summary(mod_gomp)
check_o(mod_gomp, plot = FALSE)

## ----loglik-compare-----------------------------------------------------------
mod_uni_w <- onls(density ~ Asym / (1 + exp((xmid - log(conc)) / scal)),
                  data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1),
                  sigma_x = 0.05, sigma_y = 0.1)
AIC(logLik_o(mod_uni))
AIC(logLik_o(mod_uni_w))

## ----confint------------------------------------------------------------------
set.seed(123)
confint(mod_uni, k = 100)

