## onls demo: gauge onls against published and independent reference values
## demo("gauge", package = "onls")  or  doGauge()
##
## Every comparison prints one line "item  reference  onls  |diff|  tol  RESULT", with
##   PASS  agreement within the tolerance,
##   NOTE  a documented discrepancy of the reference (explained in a comment),
##   DIFF  an undocumented deviation (expected: none).
## Kinds of reference (column 2):
##   published   value printed in the cited publication; the tolerance is half a unit of its
##               last printed decimal (the value is given as a character string); for results of
##               iterative ODRPACK runs, the convergence of that run is added (1e-3 of the standard
##               error for parameters, a relative 1e-4 for standard errors and 1e-7 for the objective;
##               helper gpub())
##   closed form exact solution computed here (SVD); relative tolerance 1e-6
##   ODRPACK     ODRPACK 2.01 (Boggs et al. 1989) through scipy.odr 1.18.1, same data, model,
##               weights and starting values, run to convergence with sstol = partol = 1e-14
##               (computed 29 October 2026); parameters are compared to 1e-3 of their standard
##               error, standard errors to a relative 1e-4 and the objective to a relative 1e-7
##   NIST        certified values of the NIST StRD nonlinear regression data sets (ordinary least
##               squares), reproduced by onls in the least-squares limit (negligible sigma_x);
##               relative tolerance 1e-6 for parameters, 1e-9 for the residual sum of squares
##   consistency properties that must hold (section 12: functions around the fit, agreement of the
##               two algorithms); PASS if the property holds
library(onls)

.GAUGE <- NULL
gauge <- function(item, value, ref, kind = "published", tol = NULL, note = FALSE, extra = 0) {
  if (is.character(ref)) {
    d <- if (grepl(".", ref, fixed = TRUE)) nchar(sub("^[^.]*\\.", "", sub("[eE].*$", "", ref))) else 0
    ex <- if (grepl("[eE]", ref)) as.numeric(sub("^.*[eE]", "", ref)) else 0
    refn <- as.numeric(ref); tol <- 0.5 * 10^(ex - d) + extra; refs <- ref
  } else { refn <- ref; refs <- formatC(ref, digits = 8, format = "g") }
  value <- unname(value); dif <- abs(value - refn); ok <- is.finite(dif) && dif <= tol
  res <- if (ok) "PASS" else if (note) "NOTE" else "DIFF"
  cat(sprintf("  %-34s %-11s %14s  onls %16s  |diff| %8.2g  tol %8.2g  %s\n",
              item, kind, refs, formatC(value, digits = 10, format = "g"), dif, tol, res))
  .GAUGE <<- rbind(.GAUGE, data.frame(item = item, kind = kind, reference = refn, onls = value,
                                      diff = dif, tol = tol, result = res))
  invisible(ok)
}
## published results of an ODRPACK run: parameter (with its SE), SE, objective
gpub <- function(item, value, ref, se = NULL, what = c("par", "se", "obj"), note = FALSE) {
  what <- match.arg(what); r <- abs(as.numeric(ref))
  gauge(item, value, ref, note = note,
        extra = switch(what, par = 1e-3 * se, se = 1e-4 * r, obj = 1e-7 * r))
}
## parameters, standard errors and objective against an ODRPACK run
gauge_odr <- function(label, m, beta, se, ss) {
  for (j in seq_along(beta))
    gauge(paste(label, names(coef(m))[j]), coef(m)[j], beta[j], "ODRPACK", tol = 1e-3 * se[j])
  for (j in seq_along(se)) if (se[j] > 0)
    gauge(paste(label, "SE", names(coef(m))[j]), m$std_errors[j], se[j], "ODRPACK", tol = 1e-4 * se[j])
  gauge(paste(label, "objective S"), m$objective, ss, "ODRPACK", tol = 1e-7 * ss)
}

######## 1. ODRPACK User's Reference Guide 2.01 (Boggs et al. 1992), example 2.C.i, pp. 39ff ########
x <- c(0, 0, 5, 7, 7.5, 10, 16, 26, 30, 34, 34.5, 100)
y <- c(1265, 1263.6, 1258, 1254, 1253, 1249.8, 1237, 1218, 1220.6, 1213.8, 1215.5, 1212)
mod1 <- onls(y ~ b1 + b2 * (exp(b3 * x) - 1)^2, data = data.frame(x, y),
             start = list(b1 = 1500, b2 = -50, b3 = -0.1))
print(mod1)
cat("Comparison with the ODRPACK guide, pp. 47-48:\n")
se1 <- mod1$std_errors
gpub("objective S", mod1$objective, "21.445", what = "obj")
gpub("b1", coef(mod1)[1], "1264.65481", se1[1]); gpub("SE b1", se1[1], "1.03492", what = "se")
gpub("b2", coef(mod1)[2], "-54.01838", se1[2]);  gpub("SE b2", se1[2], "1.583992", what = "se")
gpub("b3", coef(mod1)[3], "-0.08785", se1[3]);   gpub("SE b3", se1[3], "6.33222E-3", what = "se")

######## 2. Algorithm 676: ODRPACK (Boggs et al. 1989), pp. 355-356 and 363 ########
x <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 85, 90, 95, 100, 105)
y <- c(4.14, 8.52, 16.31, 32.18, 64.62, 98.76, 151.13, 224.74, 341.35, 423.36, 522.78, 674.32, 782.04, 920.01)
mod2 <- onls(y ~ b1 * 10^(b2 * x/(b3 + x)), data = data.frame(x, y), start = list(b1 = 1, b2 = 5, b3 = 100))
print(mod2)
cat("Comparison with Algorithm 676, p. 363:\n")
se2 <- mod2$std_errors
gpub("objective S", mod2$objective, "15.263", what = "obj")
gpub("b1", coef(mod2)[1], "4.4879", se2[1]); gpub("b2", coef(mod2)[2], "7.1882", se2[2]); gpub("b3", coef(mod2)[3], "221.8383", se2[3])
gpub("SE b1", se2[1], "0.56876", what = "se"); gpub("SE b2", se2[2], "0.69504", what = "se"); gpub("SE b3", se2[3], "37.2313", what = "se")
cat("Comparison with ODRPACK (scipy.odr), run to convergence:\n")
gauge_odr("Alg. 676", mod2, c(4.487870267, 7.188154741, 221.8377866), c(0.5687644687, 0.6950595121, 37.23230152), 15.2628142759)

######## 3. Bounded example: simple_example.f90 of ODRPACK95 (TOMS 869, Zwolak et al. 2007) ########
x <- c(0.982, 1.998, 4.978, 6.01); y <- c(2.7, 7.4, 148.0, 403.0)
## parameter bounds: the structured algorithm (active set) reaches the bounded optimum; the dense
## algorithm (minpack.lm's bounded LM, the default for n <= 100) stops slightly before it
mod3 <- onls(y ~ b1 * exp(b2 * x), data = data.frame(x, y), start = list(b1 = 2, b2 = 0.5),
             lower = c(0, 0), upper = c(10, 0.9), algorithm = "structured")
print(mod3)
cat("Comparison with ODRPACK95 and with ODRPACK (b2 held at the bound 0.9):\n")
gauge("b2 at the upper bound", coef(mod3)[2], "0.9")
## NOTE: the ODRPACK95 reference run reports b1 = 1.6334 with S = 0.2674; with b2 at its bound
## the minimum over b1 is lower, S = 0.19187 at b1 = 1.43998 (ODRPACK with b2 fixed, below)
gauge("b1 (ODRPACK95 run)", coef(mod3)[1], "1.6334", note = TRUE)
gauge("objective S (ODRPACK95 run)", mod3$objective, "0.2674", note = TRUE)
gauge("b1 (ODRPACK, b2 fixed at 0.9)", coef(mod3)[1], 1.439981541, "ODRPACK", tol = 1e-3 * 0.1672152808)
gauge("objective S (ODRPACK, b2 = 0.9)", mod3$objective, 0.191868103185, "ODRPACK", tol = 1e-7 * 0.191868103185)

######## 4. Pearson's data with York's weights (York 1966) ########
## reference values of the York et al. (2004) equations, as reported for these data by the R
## package bfsl (best-fit straight line)
x <- c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4)
y <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
sd_x <- 1/sqrt(c(1000.0, 1000.0, 500.0, 800.0, 200.0, 80.0, 60.0, 20.0, 1.8, 1.0))
sd_y <- 1/sqrt(c(1.0, 1.8, 4.0, 8.0, 20.0, 20.0, 70.0, 70.0, 100.0, 500.0))
mod4 <- onls(y ~ b0 + b1 * x, data = data.frame(x, y), start = list(b0 = 5, b1 = -0.5), sigma_x = sd_x, sigma_y = sd_y)
print(summary(mod4))
cat("Comparison with York et al. (2004), as reported by bfsl:\n")
gauge("intercept", coef(mod4)[1], "5.47991"); gauge("SE intercept", mod4$std_errors[1], "0.29497")
gauge("slope", coef(mod4)[2], "-0.48053");    gauge("SE slope", mod4$std_errors[2], "0.05799")

######## 5. Daeron & Vermeesch (2024), Table 3 / Figure 2C ########
mod5 <- onls(y ~ a + b * x, data = data.frame(x = c(9, 19, 31, 41), y = c(21, 31, 39, 49)),
             start = list(a = 10, b = 1), sigma_x = 1, sigma_y = 1)
print(mod5)
cat("Comparison with Daeron & Vermeesch (2024), Table 3:\n")
gauge("intercept a", coef(mod5)[1], "13.71")
## NOTE: with sigma_x = sigma_y = 1 this is Deming regression with variance ratio 1, whose exact
## slope is 0.85163 (closed form, below); the published 0.851 is truncated, not rounded
gauge("slope b", coef(mod5)[2], "0.851", note = TRUE)
dm <- function(x, y) { sxx <- sum((x - mean(x))^2); syy <- sum((y - mean(y))^2); sxy <- sum((x - mean(x)) * (y - mean(y)))
  b <- (syy - sxx + sqrt((syy - sxx)^2 + 4 * sxy^2)) / (2 * sxy); c(a = mean(y) - b * mean(x), b = b) }
DM5 <- dm(c(9, 19, 31, 41), c(21, 31, 39, 49))
gauge("intercept a", coef(mod5)[1], DM5[[1]], "closed form", tol = 1e-6 * abs(DM5[[1]]))
gauge("slope b", coef(mod5)[2], DM5[[2]], "closed form", tol = 1e-6 * abs(DM5[[2]]))

######## 6. Deming regression, XLSTAT help article 'Run a Deming regression to compare methods' ########
x <- c(9.8, 9.7, 10.7, 10.9, 12.4, 12.5, 12.8, 12.8, 12.9, 13.3, 13.4, 13.5, 13.7, 14.9, 15.2, 15.5)
y <- c(10.1, 11.4, 10.8, 11.3, 11.8, 12.1, 12.3, 13.6, 14.2, 14.4, 14.6, 15.3, 15.5, 15.8, 16.2, 16.5)
mod6 <- onls(y ~ a + b * x, data = data.frame(x, y), start = list(a = 2, b = 3),
             control = list(ftol = 1e-13, ptol = 1e-13))
print(mod6)
cat("Comparison with XLSTAT:\n")
gauge("intercept", coef(mod6)[1], "-1.909"); gauge("slope", coef(mod6)[2], "1.208")
DM6 <- dm(x, y)
gauge("intercept", coef(mod6)[1], DM6[[1]], "closed form", tol = 1e-6 * abs(DM6[[1]]))
gauge("slope", coef(mod6)[2], DM6[[2]], "closed form", tol = 1e-6 * abs(DM6[[2]]))

######## 7. Total least squares, closed form by SVD (Golub & Van Loan 1980) ########
set.seed(11); n <- 40
x1t <- runif(n, 0, 10); x2t <- runif(n, 0, 10)
DAT <- data.frame(x1 = x1t + rnorm(n, 0, 0.5), x2 = x2t + rnorm(n, 0, 0.5),
                  y = 3 + 1.5 * x1t - 0.8 * x2t + rnorm(n, 0, 0.5))
mod7 <- onls(y ~ b0 + b1 * x1 + b2 * x2, data = DAT, start = list(b0 = 1, b1 = 1, b2 = 1))
Z <- scale(cbind(DAT$x1, DAT$x2, DAT$y), scale = FALSE); v <- svd(Z)$v[, 3]
slope <- -v[1:2]/v[3]; TLS <- c(mean(DAT$y) - sum(slope * colMeans(DAT[, 1:2])), slope)
print(rbind(TLS = TLS, onls = coef(mod7)))
cat("Comparison with the closed-form TLS solution:\n")
for (j in 1:3) gauge(paste("TLS", names(coef(mod7))[j]), coef(mod7)[j], TLS[j], "closed form", tol = 1e-6 * abs(TLS[j]))

######## 8. Correlated predictor errors: generalized TLS, closed form by whitening + SVD ########
set.seed(123); n <- 40
Sigma <- matrix(c(0.25, 0.15, 0.15, 0.16), 2); sigma_y <- 0.3
xt <- cbind(runif(n, 0, 10), runif(n, 0, 10)); E <- matrix(rnorm(2 * n), n) %*% chol(Sigma)
DAT <- data.frame(x1 = xt[, 1] + E[, 1], x2 = xt[, 2] + E[, 2],
                  y = 3 + 1.5 * xt[, 1] - 0.8 * xt[, 2] + rnorm(n, 0, sigma_y))
mod8 <- onls(y ~ b0 + b1 * x1 + b2 * x2, data = DAT, start = list(b0 = 1, b1 = 1, b2 = 1),
             sigma_x = Sigma, sigma_y = sigma_y, control = list(ftol = 1e-13, ptol = 1e-13))
L <- chol(solve(Sigma)); U <- as.matrix(DAT[, c("x1", "x2")]) %*% t(L); w <- DAT$y/sigma_y
V <- svd(cbind(scale(U, scale = FALSE), w - mean(w)))$v[, 3]; ww <- -V[1:2]/V[3]
gTLS <- c(sigma_y * (mean(w) - sum(ww * colMeans(U))), sigma_y * drop(t(L) %*% ww))
print(rbind(gTLS = gTLS, onls = coef(mod8)))
cat("Comparison with the closed-form generalized TLS solution:\n")
for (j in 1:3) gauge(paste("gTLS", names(coef(mod8))[j]), coef(mod8)[j], gTLS[j], "closed form", tol = 1e-6 * abs(gTLS[j]))

######## 9. ODRPACK weights WE (response) and WD (predictor), per observation ########
set.seed(123); n <- 30
xt <- seq(0.5, 10, length.out = n); WD <- runif(n, 0.5, 4); WE <- runif(n, 0.5, 4)
x <- xt + rnorm(n, 0, 0.3/sqrt(WD)); y <- 2 * exp(-0.3 * xt) + 0.5 + rnorm(n, 0, 0.03/sqrt(WE))
mod9 <- onls(y ~ a * exp(-b * x) + c, data = data.frame(x, y), start = list(a = 1.5, b = 0.2, c = 0.3),
             weights = WE, sigma_x = 1/sqrt(WD), known_sigma = FALSE)
print(summary(mod9))
cat("Comparison with ODRPACK (we = WE, wd = WD):\n")
gauge_odr("WE/WD", mod9, c(2.046024826, 0.3048444214, 0.5055123772), c(0.02444419649, 0.01170446794, 0.02096574718), 0.0521643094555)

######## 10. NIST StRD problems as orthogonal regressions (unit precisions) ########
## NIST certifies these problems for ordinary least squares only; the orthogonal fits are gauged
## against ODRPACK with the same starting values (ENSO: from the certified values, because its
## periodic terms have many local minima and ODRPACK and onls end in different ones from the
## NIST starting values).
## Chwirut2
x <- c(0.5, 1, 1.75, 3.75, 5.75, 0.875, 2.25, 3.25, 5.25, 0.75, 1.75, 2.75,
       4.75, 0.625, 1.25, 2.25, 4.25, 0.5, 3, 0.75, 3, 1.5, 6, 3,
       6, 1.5, 3, 0.5, 2, 4, 0.75, 2, 5, 0.75, 2.25, 3.75,
       5.75, 3, 0.75, 2.5, 4, 0.75, 2.5, 4, 0.75, 2.5, 4, 0.5,
       6, 3, 0.5, 2.75, 0.5, 1.75)
y <- c(92.9, 57.1, 31.05, 11.5875, 8.025, 63.6, 21.4, 14.25, 8.475, 63.8,
       26.8, 16.4625, 7.125, 67.3, 41, 21.15, 8.175, 81.5, 13.12, 59.9,
       14.62, 32.9, 5.44, 12.56, 5.44, 32, 13.95, 75.8, 20, 10.42,
       59.5, 21.67, 8.55, 62, 20.2, 7.76, 3.75, 11.81, 54.7, 23.7,
       11.55, 61.3, 17.7, 8.74, 59.2, 16.3, 8.62, 81, 4.87, 14.62,
       81.7, 17.17, 81.3, 28.9)
CHW <- data.frame(x, y)
m_chw <- onls(y ~ exp(-b1 * x)/(b2 + b3 * x), data = CHW, start = list(b1 = 0.1, b2 = 0.005, b3 = 0.01))
print(m_chw); plot(m_chw, main = "Chwirut2")
gauge_odr("Chwirut2", m_chw, c(0.1540697415, 0.004996924934, 0.01261378928), c(0.02830538839, 0.002014073345, 0.002064198982), 5.57234046344)
## Gauss1 (full-precision NIST data; the copy in NIST() of onls <= 0.2-1 was rounded to 4 decimals)
GAU <- data.frame(x = 1:250, y = c(97.62227, 97.80724, 96.62247, 92.59022, 91.23869, 95.32704, 90.3504, 89.46235, 91.7252, 89.86916,
       86.88076, 85.9436, 87.60686, 86.25839, 80.74976, 83.03551, 88.25837, 82.01316, 82.74098, 83.30034,
       81.2785, 81.85506, 80.75195, 80.09573, 81.07633, 78.81542, 78.38596, 79.93386, 79.48474, 79.95942,
       76.10691, 78.3983, 81.4306, 82.48867, 81.65462, 80.84323, 88.68663, 84.74438, 86.83934, 85.97739,
       91.28509, 97.22411, 93.51733, 94.10159, 101.9176, 98.43134, 110.4214, 107.6628, 111.7288, 116.5115,
       120.7609, 123.9553, 124.2437, 130.7996, 133.296, 130.7788, 132.0565, 138.6584, 142.9252, 142.7215,
       144.1249, 147.4377, 148.2647, 152.0519, 147.3863, 149.2074, 148.9537, 144.5876, 148.1226, 148.0144,
       143.8893, 140.9088, 143.4434, 139.3938, 135.9878, 136.3927, 126.7262, 124.4487, 122.8647, 113.8557,
       113.7037, 106.8407, 107.0034, 102.4629, 96.09296, 94.57555, 86.98824, 84.90154, 81.18023, 76.40117,
       67.092, 72.67155, 68.10848, 67.99088, 63.34094, 60.55253, 56.18687, 53.64482, 53.70307, 48.07893,
       42.21258, 45.65181, 41.69728, 41.24946, 39.21349, 37.71696, 36.68395, 37.30393, 37.43277, 37.45012,
       32.64648, 31.84347, 31.39951, 26.68912, 32.25323, 27.61008, 33.58649, 28.10714, 30.26428, 28.01648,
       29.11021, 23.02099, 25.65091, 28.50295, 25.23701, 26.13828, 33.5326, 29.25195, 27.09847, 26.52999,
       25.52401, 26.69218, 24.55269, 27.71763, 25.20297, 25.61483, 25.06893, 27.6393, 24.94851, 25.86806,
       22.48183, 26.90045, 25.39919, 17.90614, 23.76039, 25.89689, 27.64231, 22.86101, 26.47003, 23.72888,
       27.54334, 30.52683, 28.07261, 34.92815, 28.29194, 34.19161, 35.41207, 37.09336, 40.9833, 39.53923,
       47.80123, 47.46305, 51.04166, 54.58065, 57.53001, 61.42089, 62.79032, 68.51455, 70.23053, 74.42776,
       76.59911, 81.62053, 83.42208, 79.17451, 88.56985, 85.66525, 86.55502, 90.65907, 84.2729, 85.7222,
       83.10702, 82.16884, 80.42568, 78.15692, 79.79691, 77.84378, 74.50327, 71.57289, 65.88031, 65.01385,
       60.19582, 59.66726, 52.95478, 53.87792, 44.91274, 41.09909, 41.68018, 34.53379, 34.86419, 33.14787,
       29.58864, 27.29462, 21.91439, 19.08159, 24.9029, 19.82341, 16.75551, 18.24558, 17.23549, 16.34934,
       13.71285, 14.75676, 13.97169, 12.42867, 14.35519, 7.703309, 10.23441, 11.78315, 13.87768, 4.5357,
       10.05928, 8.424824, 10.53312, 9.602255, 7.877514, 6.258121, 8.899865, 7.877754, 12.51191, 10.66205,
       6.0354, 6.790655, 8.783535, 4.600288, 8.400915, 7.216561, 10.01741, 7.331278, 6.527863, 2.842001,
       10.32507, 4.790995, 8.377101, 6.264445, 2.706213, 8.362329, 8.983658, 3.362571, 1.182746, 4.875359))
fG <- y ~ b1 * exp(-b2 * x) + b3 * exp(-(x - b4)^2/b5^2) + b6 * exp(-(x - b7)^2/b8^2)
m_gau <- onls(fG, data = GAU, start = list(b1 = 100, b2 = 0.01, b3 = 100, b4 = 50, b5 = 20, b6 = 70, b7 = 200, b8 = 20))
print(m_gau); plot(m_gau, main = "Gauss1")
gauge_odr("Gauss1", m_gau, c(99.06697747, 0.01054069446, 100.7187587, 67.52904935, 23.00483836, 72.42609706, 178.8537009, 18.40983674), c(0.5631277489, 9.757817266e-05, 0.7109993175, 0.2165940227, 0.2543925702, 0.7711102092, 0.2408370682, 0.2841064898), 714.943905161)
## ENSO
ENS <- data.frame(x = 1:168, y = c(12.9, 11.3, 10.6, 11.2, 10.9, 7.5, 7.7, 11.7, 12.9, 14.3, 10.9, 13.7, 17.1, 14,
       15.3, 8.5, 5.7, 5.5, 7.6, 8.6, 7.3, 7.6, 12.7, 11, 12.7, 12.9, 13, 10.9,
       10.4, 10.2, 8, 10.9, 13.6, 10.5, 9.2, 12.4, 12.7, 13.3, 10.1, 7.8, 4.8, 3,
       2.5, 6.3, 9.7, 11.6, 8.6, 12.4, 10.5, 13.3, 10.4, 8.1, 3.7, 10.7, 5.1, 10.4,
       10.9, 11.7, 11.4, 13.7, 14.1, 14, 12.5, 6.3, 9.6, 11.7, 5, 10.8, 12.7, 10.8,
       11.8, 12.6, 15.7, 12.6, 14.8, 7.8, 7.1, 11.2, 8.1, 6.4, 5.2, 12, 10.2, 12.7,
       10.2, 14.7, 12.2, 7.1, 5.7, 6.7, 3.9, 8.5, 8.3, 10.8, 16.7, 12.6, 12.5, 12.5,
       9.8, 7.2, 4.1, 10.6, 10.1, 10.1, 11.9, 13.6, 16.3, 17.6, 15.5, 16, 15.2, 11.2,
       14.3, 14.5, 8.5, 12, 12.7, 11.3, 14.5, 15.1, 10.4, 11.5, 13.4, 7.5, 0.6, 0.3,
       5.5, 5, 4.6, 8.2, 9.9, 9.2, 12.5, 10.9, 9.9, 8.9, 7.6, 9.5, 8.4, 10.7,
       13.6, 13.7, 13.7, 16.5, 16.8, 17.1, 15.4, 9.5, 6.1, 10.1, 9.3, 5.3, 11.2, 16.6,
       15.6, 12, 11.5, 8.6, 13.8, 8.7, 8.6, 8.6, 8.7, 12.8, 13.2, 14, 13.4, 14.8))
fE <- y ~ b1 + b2 * cos(2 * pi * x/12) + b3 * sin(2 * pi * x/12) + b5 * cos(2 * pi * x/b4) +
  b6 * sin(2 * pi * x/b4) + b8 * cos(2 * pi * x/b7) + b9 * sin(2 * pi * x/b7)
certE <- c(b1 = 1.0510749193E+01, b2 = 3.0762128085E+00, b3 = 5.3280138227E-01, b4 = 4.4311088700E+01,
           b5 = -1.6231428586E+00, b6 = 5.2554493756E-01, b7 = 2.6887614440E+01, b8 = 2.1232288488E-01,
           b9 = 1.4966870418E+00)
m_ens <- onls(fE, data = ENS, start = as.list(certE))
print(m_ens); plot(m_ens, main = "ENSO")
gauge_odr("ENSO", m_ens, c(10.35444359, 5.668019513, 0.998402262, 43.63303044, -1.922521051, 0.2637078762, 26.65269076, 0.04440235759, 1.723689926), c(0.2362854023, 0.3220735025, 0.3187907635, 1.116417333, 0.3337904403, 0.6905058541, 0.4761326501, 0.6708714995, 0.3282527109), 244.133615034)
## Rat43 data with the 3-parameter Gompertz curve (the 4-parameter Rat43 model degenerates to it
## in the orthogonal fit)
RAT <- data.frame(x = 1:15, y = c(16.08, 33.83, 65.8, 97.2, 191.55, 326.2, 386.87, 520.53,
       590.03, 651.92, 724.93, 699.56, 689.96, 637.56, 717.41))
m_rat <- onls(y ~ b1 * exp(-exp(c - b3 * x)), data = RAT, start = list(b1 = 800, c = 1, b3 = 0.2))
print(m_rat); plot(m_rat, main = "Rat43 data, Gompertz model")
gauge_odr("Rat43 Gompertz", m_rat, c(837.3611933, 1.716313049, 0.2755950805), c(101.9243997, 0.2782533807, 0.06801957186), 23.421333659)

######## 11. NIST certified values in the least-squares limit (sigma_x -> 0) ########
## With a negligible predictor error, the orthogonal fit becomes ordinary least squares, so onls
## must reproduce the NIST certified parameters and residual sums of squares.
ols <- function(form, d, start, sx)
  onls(form, data = d, start = start, sigma_x = sx, known_sigma = FALSE,
       control = list(ftol = 1e-15, ptol = 1e-15))
nist <- function(label, m, cert, rss) {
  for (j in seq_along(cert)) gauge(paste(label, names(cert)[j]), coef(m)[j], cert[j], "NIST", tol = 1e-6 * abs(cert[j]))
  gauge(paste(label, "RSS"), deviance(m), rss, "NIST", tol = 1e-9 * rss)
}
certC <- c(b1 = 1.6657666537E-01, b2 = 5.1653291286E-03, b3 = 1.2150007096E-02)
nist("Chwirut2 (OLS)", ols(y ~ exp(-b1 * x)/(b2 + b3 * x), CHW, as.list(certC * 0.9), 1e-7), certC, 5.1304802941E+02)
certG <- c(b1 = 9.8778210871E+01, b2 = 1.0497276517E-02, b3 = 1.0048990633E+02, b4 = 6.7481111276E+01,
           b5 = 2.3129773360E+01, b6 = 7.1994503004E+01, b7 = 1.7899805021E+02, b8 = 1.8389389025E+01)
nist("Gauss1 (OLS)", ols(fG, GAU, list(b1 = 97, b2 = 0.009, b3 = 100, b4 = 65, b5 = 20, b6 = 70, b7 = 178, b8 = 16.5), 1e-6),
     certG, 1.3158222432E+03)
nist("ENSO (OLS)", ols(fE, ENS, as.list(certE * c(0.99, 1.01, 0.99, 1.001, 1.01, 0.99, 1.001, 0.99, 1.01)), 1e-6),
     certE, 7.8853978668E+02)
certR <- c(b1 = 6.9964151270E+02, b2 = 5.2771253025E+00, b3 = 7.5962938329E-01, b4 = 1.2792483859E+00)
nist("Rat43 (OLS)", ols(y ~ b1/((1 + exp(b2 - b3 * x))^(1/b4)), RAT, list(b1 = 700, b2 = 5, b3 = 0.75, b4 = 1.3), 1e-6),
     certR, 8.7864049080E+03)

######## 12. Consistency checks of the functions around the fit ########
## Not a comparison with external values: each line checks a property that must hold (bugs fixed
## in 0.2-1, and the agreement of the two algorithms). PASS if it holds, DIFF otherwise.
gcheck <- function(item, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) FALSE)
  res <- if (ok) "PASS" else "DIFF"
  cat(sprintf("  %-72s %-11s  %s\n", item, "consistency", res))
  .GAUGE <<- rbind(.GAUGE, data.frame(item = item, kind = "consistency", reference = 1, onls = as.numeric(ok),
                                      diff = as.numeric(!ok), tol = 0, result = res))
  invisible(ok)
}
## warning (with a pattern) or error raised?
warns <- function(expr, pattern) {
  w <- character(0)
  withCallingHandlers(tryCatch(expr, error = function(e) NULL),
                      warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  any(grepl(pattern, w))
}
fails <- function(expr, pattern) tryCatch({ expr; FALSE }, error = function(e) grepl(pattern, conditionMessage(e)))
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
## both algorithms from the same call: same parameters, objective and SEs, converged, orthogonal
same_fit <- function(label, call, tol = 1e-6) {
  cs <- call; cs$algorithm <- "structured"; cd <- call; cd$algorithm <- "dense"
  ms <- quiet(eval(cs, parent.frame())); md <- quiet(eval(cd, parent.frame()))
  rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))
  gcheck(paste(label, ": structured = dense, parameters"), rel(coef(ms), coef(md)) < tol)
  gcheck(paste(label, ": structured = dense, objective"), rel(ms$objective, md$objective) < tol)
  gcheck(paste(label, ": structured = dense, SEs"), rel(ms$std_errors[ms$std_errors > 0], md$std_errors[md$std_errors > 0]) < 10 * tol)
  gcheck(paste(label, ": structured converged, all points orthogonal"), ms$convergence == 0 && all(ms$ortho))
}

DN <- subset(DNase, Run == 1); set.seed(123)
DN$density <- sapply(DN$density, function(x) rnorm(1, x, 0.1 * x))
fDN <- density ~ Asym / (1 + exp((xmid - log(conc)) / scal)); stDN <- list(Asym = 3, xmid = 0, scal = 1)
mDN <- onls(fDN, DN, start = stDN)
xPY <- c(0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4); yPY <- c(5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5)
sxPY <- 1/sqrt(c(1000, 1000, 500, 800, 200, 80, 60, 20, 1.8, 1.0)); syPY <- 1/sqrt(c(1, 1.8, 4, 8, 20, 20, 70, 70, 100, 500))
DPY <- data.frame(x = xPY, y = yPY)
set.seed(5); n <- 40; S2 <- matrix(c(0.04, 0.024, 0.024, 0.09), 2)
xt <- cbind(runif(n, 1, 5), runif(n, 1, 5)); E <- matrix(rnorm(2 * n), n) %*% chol(S2)
DC <- data.frame(x1 = xt[, 1] + E[, 1], x2 = xt[, 2] + E[, 2], y = 2 * exp(0.3 * xt[, 1]) + 1.5 * log(xt[, 2]) + rnorm(n, 0, 0.2))
fDC <- y ~ b1 * exp(b2 * x1) + b3 * log(x2); stDC <- list(b1 = 1, b2 = 0.2, b3 = 1)

cat("Inference and helper functions:\n")
m1 <- mDN; m1$ortho[1] <- FALSE; set.seed(1)
ci <- quiet(confint(m1, method = "boot", k = 5, max_tries = 20))
gcheck("confint(boot): stops after max_tries, k fits collected", attr(ci, "n_boot") == 5)
mPY <- onls(y ~ b0 + b1 * x, DPY, start = list(b0 = 5, b1 = -0.5), sigma_x = sxPY, sigma_y = syPY)
idx <- c(5, 6, 6, 8, 1, 1, 9, 2, 1, 3)
man <- onls(y ~ b0 + b1 * x, DPY[idx, ], start = list(b0 = 5, b1 = -0.5), sigma_x = sxPY[idx], sigma_y = syPY[idx])
gcheck("bootstrap refit resamples sigma_x and sigma_y with their rows", isTRUE(all.equal(coef(onls:::.refit(mPY, idx = idx)), coef(man), tolerance = 1e-6)))
mDC <- onls(fDC, DC, start = stDC, sigma_x = S2, sigma_y = 0.2)
gcheck("check_o: correlated sigma_x fit, all points orthogonal", all(check_o(mDC, plot = FALSE)$Ortho))
set.seed(1); DNs <- DN[sample(nrow(DN)), ]; DNs$density[c(2, 5)] <- NA
mNA <- onls(fDN, DNs, start = stDN, na.action = na.exclude); okNA <- !is.na(DNs$density)
gcheck("x0()/y0() in data order, NA re-inserted with na.exclude",
       identical(which(is.na(x0(mNA))), c(2L, 5L)) && isTRUE(all.equal(unname(x0(mNA)[okNA]), unname(mNA$xi[, 1]))) &&
       isTRUE(all.equal(unname(y0(mNA)[okNA]), unname(mNA$yi))) && all(abs(x0(mNA)[okNA] - DNs$conc[okNA]) < 0.5))
gcheck("check_o in data order; all points orthogonal with na.exclude", isTRUE(all.equal(check_o(mNA, plot = FALSE)$conc, DNs$conc[okNA])) && all(mNA$ortho))
mFX <- onls(fDN, DN, start = stDN, fixed = c(TRUE, FALSE, FALSE))
gcheck("predict() with a fixed parameter: confidence interval", identical(dim(predict(mFX, newdata = data.frame(conc = c(1, 6)), interval = "confidence")), c(2L, 3L)))
gcheck("summary(correlation = TRUE) with a fixed parameter", identical(dim(summary(mFX, correlation = TRUE)$correlation), c(2L, 2L)))
pr <- profile(mDN); cip <- quiet(confint(pr))
gcheck("profile() dispatches to profile.onls", inherits(pr, "profile.onls"))
gcheck("confint(profile): finite limits around the estimates", !anyNA(cip) && all(cip[, 1] < coef(mDN) & coef(mDN) < cip[, 2]))
gcheck("confint(fit) defaults to the profile interval", isTRUE(all.equal(quiet(confint(mDN)), cip, tolerance = 0.01)))
gcheck("plot(profile) runs", { plot(pr); TRUE })
g <- function() { mylog <- function(z) log(z); onls(density ~ Asym / (1 + exp((xmid - mylog(conc)) / scal)), DN, start = stDN) }
mG <- g()
gcheck("predict() finds functions of the formula environment", length(predict(mG, newdata = data.frame(conc = 6))) == 1)
gcheck("plot() finds functions of the formula environment", { plot(mG); TRUE })
h <- function() { D <- DN; onls(fDN, D, start = stDN) }
set.seed(1)
gcheck("confint(boot) for a fit made inside a function", identical(dim(quiet(confint(h(), method = "boot", k = 5))), c(3L, 2L)))
mNO <- mDN; mNO$ortho[] <- NA
gcheck("print() with a missing orthogonality check", any(grepl("could not be checked", capture.output(print(mNO)))))
set.seed(1); xw <- 1:40; yw <- 3 + 2 * xw + rnorm(40, 0, 3)
gcheck("warning when extend/window foot points stop early",
       warns(onls(y ~ a + b * x, data.frame(x = xw, y = yw), start = list(a = 1, b = 1), extend = c(0.1, 0.1), window = 3), "bounded solver"))
gcheck("warning when sigma_y is given without sigma_x", warns(onls(fDN, DN, start = stDN, sigma_y = 0.1), "sigma_x"))
set.seed(3); n <- 30; tt <- 1:n; z <- runif(n, 1, 5); yt <- 1 + 0.2 * tt + 2 * sqrt(z) + rnorm(n, 0, 0.1); z <- z + rnorm(n, 0, 0.1)
gcheck("tiny sigma_x (error-free predictor): no false orthogonality failures",
       all(onls(y ~ a + b * tt + c * sqrt(z), data.frame(tt, z, y = yt), start = list(a = 1, b = 1, c = 1), sigma_x = c(1e-8, 0.1), sigma_y = 0.1)$ortho))
wz <- rep(1, nrow(DN)); wz[3] <- 0; mW <- onls(fDN, DN, start = stDN, weights = wz)
gcheck("zero weight: logLik_o finite, df_resid without that observation", is.finite(logLik_o(mW)) && mW$df_resid == nrow(DN) - 1L - 3L)
set.seed(2); xz <- c(0, 0, 1:18); yz <- 2 + 0.5 * xz + 0.02 * xz^2 + rnorm(20, 0, 0.3); wz2 <- rep(1, 20); wz2[1] <- 0
gcheck("zero-weight observation at x = 0 counts as orthogonal",
       all(onls(y ~ a + b * x + c * x^2, data.frame(x = xz, y = yz), start = list(a = 1, b = 1, c = 0), weights = wz2)$ortho))
gcheck("NA in sigma_y gives an error",
       fails(onls(y ~ a + b * x, data.frame(x = 1:5, y = c(1, 2, 3.1, 4, 5.2)), start = list(a = 0, b = 1), sigma_y = c(0.1, NA, 0.1, 0.1, 0.1)), "NA"))
gcheck("confint(parm = 1:2) selects parameters by index", identical(rownames(confint(mDN, parm = 1:2, method = "wald")), c("Asym", "xmid")))
set.seed(1); x1p <- 1:15; y1p <- 3 * exp(0.1 * x1p) + rnorm(15, 0, 0.2)
m1p <- onls(y ~ a * exp(0.1 * x), data.frame(x = x1p, y = y1p), start = list(a = 3)); ci1 <- quiet(confint(m1p))
gcheck("profile interval for a model with one free parameter", !anyNA(ci1) && ci1[1, 1] < coef(m1p) && coef(m1p) < ci1[1, 2])
D3p <- data.frame(x = c(1, 2, 3), y = c(2.1, 3.9, 6.2)); m3p <- onls(y ~ a * x^b, D3p, start = list(a = 2, b = 1)); ci3 <- quiet(confint(m3p))
gcheck("profile limits ordered around the estimates (n = 3)", all(ci3[, 1] <= coef(m3p) & coef(m3p) <= ci3[, 2]))
gcheck("profile of an exact fit is refused", fails(profile(onls(y ~ a * x^b + c, D3p, start = list(a = 2, b = 1, c = 0))), "exact fit"))
set.seed(1); xb <- 1:20; yb <- 3 * exp(0.1 * xb) + rnorm(20, 0, 0.3)
mB <- onls(y ~ a * exp(b * x), data.frame(x = xb, y = yb), start = list(a = 3.1, b = 0.1), lower = c(3.05, 0))
ciB <- NULL; wB <- warns(ciB <- confint(mB), "bound")
gcheck("profile limit beyond a lower bound is NA, with a warning", wB && is.na(ciB["a", 1]) && !is.na(ciB["a", 2]))

cat("The two algorithms:\n")
same_fit("DNase", quote(onls(fDN, DN, start = stDN)))
same_fit("DNase, Asym fixed", quote(onls(fDN, DN, start = stDN, fixed = c(TRUE, FALSE, FALSE))))
same_fit("Pearson-York sigmas", quote(onls(y ~ b0 + b1 * x, DPY, start = list(b0 = 5, b1 = -0.5), sigma_x = sxPY, sigma_y = syPY)))
same_fit("correlated 2 x 2 sigma_x", quote(onls(fDC, DC, start = stDC, sigma_x = S2, sigma_y = 0.2)))
set.seed(6); n <- 40; sx6 <- cbind(runif(n, 0.05, 0.3), runif(n, 0.05, 0.3)); sy6 <- runif(n, 0.1, 0.4); xt6 <- cbind(runif(n, 1, 5), runif(n, 1, 5))
D6 <- data.frame(x1 = xt6[, 1] + rnorm(n, 0, sx6[, 1]), x2 = xt6[, 2] + rnorm(n, 0, sx6[, 2]), y = 3 * xt6[, 1] / (1 + 0.4 * xt6[, 2]) + rnorm(n, 0, sy6))
same_fit("per-observation sigma_x matrix", quote(onls(y ~ b1 * x1 / (1 + b2 * x2), D6, start = list(b1 = 2, b2 = 0.2), sigma_x = sx6, sigma_y = sy6)))
set.seed(11); n <- 60; S3 <- matrix(c(0.04, 0.02, 0.01, 0.02, 0.09, 0.03, 0.01, 0.03, 0.06), 3)
xt3 <- cbind(runif(n, 1, 5), runif(n, 1, 5), runif(n, 1, 5)); E3 <- matrix(rnorm(3 * n), n) %*% chol(S3)
D3 <- data.frame(x1 = xt3[, 1] + E3[, 1], x2 = xt3[, 2] + E3[, 2], x3 = xt3[, 3] + E3[, 3],
                 y = 1 + 0.5 * xt3[, 1] + log(xt3[, 2]) + 0.2 * xt3[, 3]^2 + rnorm(n, 0, 0.2))
same_fit("correlated 3 x 3 sigma_x", quote(onls(y ~ b0 + b1 * x1 + b2 * log(x2) + b3 * x3^2, D3,
                                           start = list(b0 = 1, b1 = 1, b2 = 1, b3 = 0.1), sigma_x = S3, sigma_y = 0.2)))
set.seed(7); n <- 40; xt7 <- seq(0.1, 2, length = n)
D7 <- data.frame(x = xt7 + rnorm(n, 0, 0.02), y = 1e7 * exp(1.5 * xt7) * (1 + rnorm(n, 0, 0.02)))
same_fit("y much more precise than x", quote(onls(y ~ A * exp(k * x), D7, start = list(A = 1e7, k = 1.5))))
set.seed(5); n <- 30; xt8 <- seq(0, 2, length = n)
D8 <- data.frame(x = xt8 + rnorm(n, 0, 0.01), y = 1e-9 * exp(2 * xt8) * (1 + rnorm(n, 0, 0.01)))
same_fit("parameters on very different scales", quote(onls(y ~ a * exp(b * x), D8, start = list(a = 1e-9, b = 2), sigma_x = 0.01, sigma_y = 1e-11)))
set.seed(3); n <- 50; xa <- rnorm(n); xb2 <- -0.9 * xa + sqrt(0.19) * rnorm(n); ya <- -xa - 0.1 * xb2 + rnorm(n)
DA <- data.frame(x1 = xa, x2 = xb2, y = ya)
mA <- onls(y ~ a * x1 + b * x2, DA, start = list(a = 0.5, b = 0.5), lower = c(0, 0), algorithm = "structured")
pA <- optimize(function(b) onls(y ~ a * x1 + b * x2, DA, start = list(a = 0, b = b), fixed = c(TRUE, TRUE))$objective, c(0, 5))
gcheck("structured, two active bounds: box optimum reached", mA$objective <= pA$objective * (1 + 1e-6) && all(mA$ortho))
xe6 <- c(0.982, 1.998, 4.978, 6.01); ye6 <- c(2.7, 7.4, 148.0, 403.0)
mS6 <- onls(y ~ b1 * exp(b2 * x), data.frame(x = xe6, y = ye6), start = list(b1 = 2, b2 = 0.5), lower = c(0, 0), upper = c(10, 0.9), algorithm = "structured")
gcheck("structured, one bound: objective <= dense", mS6$objective <= quiet(update(mS6, algorithm = "dense"))$objective * (1 + 1e-8))
set.seed(1); n <- 120; xt9 <- seq(0, 10, length.out = n)
D9 <- data.frame(x = xt9 + rnorm(n, 0, 0.2), y = 5 / (1 + exp(-(xt9 - 5))) + rnorm(n, 0, 0.2))
f9 <- y ~ b1 / (1 + exp(-(x - b2) / b3)); st9 <- list(b1 = 4, b2 = 4, b3 = 1.5)
gcheck("auto: structured for n = 120, dense for n = 100",
       onls(f9, D9, start = st9)$algorithm == "structured" && onls(f9, D9[1:100, ], start = st9)$algorithm == "dense")
set.seed(1); xq <- 1:20; yq <- 10 + 3 * xq^2 + rnorm(20, 0, 50)
gcheck("extend/window use the dense algorithm",
       onls(y ~ a + b * x^2, data.frame(x = xq, y = yq), start = list(a = 10, b = 3), extend = c(0.2, 0), algorithm = "structured")$algorithm == "dense")
set.seed(1); n <- 5000; xt10 <- seq(0, 10, length.out = n)
D10 <- data.frame(x = xt10 + rnorm(n, 0, 0.2), y = 5 / (1 + exp(-(xt10 - 5))) + rnorm(n, 0, 0.2))
tt10 <- system.time(m10 <- onls(f9, D10, start = st9))[["elapsed"]]
gcheck(sprintf("structured, n = 5000: converged in %.2f s, close to the true parameters", tt10),
       m10$convergence == 0 && tt10 < 20 && max(abs(coef(m10) - c(5, 5, 1)) / c(5, 5, 1)) < 0.05)

cat("Likelihood, summary, predict (0.3-0):\n")
llPY <- logLik_o(mPY); Np <- nrow(DPY)
llman <- -mPY$objective / 2 - Np * 2 / 2 * log(2 * pi) + 0.5 * sum(log(1 / syPY^2)) + 0.5 * sum(log(1 / sxPY^2))
gcheck("logLik_o, known variances: EIV formula, df = free parameters", abs(llPY - llman) < 1e-8 && attr(llPY, "df") == 2)
mU1 <- onls(fDN, DN, start = stDN, sigma_x = 0.05, sigma_y = 0.1, known_sigma = FALSE)
mU2 <- onls(fDN, DN, start = stDN, sigma_x = 0.5, sigma_y = 1, known_sigma = FALSE)
gcheck("logLik_o, estimated scale: invariant to a common sigma factor, df = q + 1",
       abs(logLik_o(mU1) - logLik_o(mU2)) < 1e-6 && attr(logLik_o(mU1), "df") == 4)
gcheck("summary: z values for known variances, t values otherwise",
       identical(colnames(coef(summary(mPY)))[3:4], c("z value", "Pr(>|z|)")) &&
       identical(colnames(coef(summary(mDN)))[3:4], c("t value", "Pr(>|t|)")))
ndPY <- data.frame(x = c(1, 4, 7)); pc <- predict(mPY, ndPY, interval = "confidence", se.fit = TRUE)
gcheck("predict, known variances: normal quantiles", isTRUE(all.equal(unname(pc$fit[, "upr"] - pc$fit[, "fit"]), qnorm(0.975) * pc$se.fit)) && is.infinite(pc$df))
pp <- predict(mPY, ndPY, interval = "prediction", sigma_x = 0.2, sigma_y = 0.3)
gcheck("predict, line: prediction se^2 = se_fit^2 + sy^2 + b1^2 sx^2",
       isTRUE(all.equal(unname(pp[, "upr"] - pp[, "fit"]), qnorm(0.975) * sqrt(pc$se.fit^2 + 0.3^2 + coef(mPY)[["b1"]]^2 * 0.2^2))))
pd <- predict(mDN, data.frame(conc = c(1, 5)), interval = "prediction", se.fit = TRUE)
gcheck("predict, estimated scale: t quantiles with df_resid", pd$df == mDN$df_resid && isTRUE(all.equal(pd$residual.scale^2, mDN$reduced_chisq)))
ndC <- DC[1:4, ]
gcheck("predict, correlated sigma_x: default = covariance given for newdata",
       isTRUE(all.equal(predict(mDC, ndC, interval = "prediction"), predict(mDC, ndC, interval = "prediction", sigma_x = S2))))
gcheck("predict: p x p sigma_x is a covariance also when nrow(newdata) = p",
       isTRUE(all.equal(predict(mDC, DC[1:2, ], interval = "prediction"), predict(mDC, DC[1:3, ], interval = "prediction")[1:2, ])))
psc <- predict(mDN, data.frame(conc = c(1, 5)), se.fit = TRUE, scale = 2 * sqrt(mDN$reduced_chisq))
gcheck("predict(scale =) rescales se.fit, as predict.lm", isTRUE(all.equal(psc$se.fit, 2 * pd$se.fit)))
gcheck("predict: per-observation sigma_x of the fit needs sigma_x for newdata",
       fails(predict(mPY, ndPY, interval = "prediction"), "sigma_x|sigma_y"))

cat("Arguments evaluated in data (nls-style), algorithm choice, representation:\n")
DNw <- DN; DNw$w <- seq(0.5, 1.5, length.out = nrow(DN))
mE <- onls(fDN, DNw, start = stDN, weights = w, subset = conc > 0.1)
keep <- DNw$conc > 0.1
mM <- onls(fDN, DNw[keep, ], start = stDN, weights = DNw$w[keep])
gcheck("weights and subset evaluated in data, as in nls()", isTRUE(all.equal(coef(mE), coef(mM), tolerance = 1e-8)) && mE$n == sum(keep))
DPYs <- data.frame(DPY, sx = sxPY, sy = syPY)
DNw$flag <- DNw$conc > 0.1; DNw$flag[3] <- NA
gcheck("NA in a logical subset drops the observation, as in nls()",
       isTRUE(all.equal(coef(onls(fDN, DNw, start = stDN, weights = w, subset = flag)),
                        coef(onls(fDN, DNw[which(DNw$flag), ], start = stDN, weights = w)), tolerance = 1e-8)))
gcheck("sigma_x and sigma_y as columns of data", isTRUE(all.equal(coef(onls(y ~ b0 + b1 * x, DPYs, start = list(b0 = 5, b1 = -0.5), sigma_x = sx, sigma_y = sy)), coef(mPY), tolerance = 1e-8)))
gcheck("auto: structured when parameter bounds are given (n <= 100)",
       onls(fDN, DN, start = stDN, lower = c(0, -Inf, 0))$algorithm == "structured")
Lq <- matrix(mDC$Lx[1, , ], 2)
gcheck("Lx: t(L) L = solve(Sigma_x) for a p x p covariance", isTRUE(all.equal(crossprod(Lq), solve(S2), check.attributes = FALSE)))
sdv <- c(0.15, 0.25); mV <- onls(fDC, DC, start = stDC, sigma_x = sdv, sigma_y = 0.2)
mR <- onls(fDC, DC, start = stDC, sigma_x = matrix(sdv, nrow(DC), 2, byrow = TRUE), sigma_y = 0.2)
mS <- onls(fDC, DC, start = stDC, sigma_x = diag(sdv^2), sigma_y = 0.2)
gcheck("same fit: sigma_x as length-p vector, n x p matrix, p x p covariance",
       isTRUE(all.equal(coef(mV), coef(mR), tolerance = 1e-8)) && isTRUE(all.equal(coef(mV), coef(mS), tolerance = 1e-8)) &&
       isTRUE(all.equal(mV$std_errors, mS$std_errors, tolerance = 1e-6)))
gcheck("fit object: 36 elements, no duplicated components of 0.2",
       length(mDN) == 36 && !any(c("parONLS", "x0", "resid_o", "pred", "QR", "Q_x", "model") %in% names(mDN)))

######## Summary ########
cat("\nSummary of all comparisons:\n")
print(table(factor(.GAUGE$result, levels = c("PASS", "NOTE", "DIFF"))))
