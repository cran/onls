### R code from vignette source 'onls.Rnw'

###################################################
### code chunk number 1: onls.Rnw:58-63
###################################################
library(onls)
DNase1 <- subset(DNase, Run == 1)
DNase1$density <- sapply(DNase1$density, function(x) rnorm(1, x, 0.1 * x))
mod1 <- onls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)), 
             data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1))


###################################################
### code chunk number 2: onls.Rnw:76-77
###################################################
mod1


###################################################
### code chunk number 3: onls.Rnw:86-92
###################################################
x <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 85, 90, 95, 100, 105)
y <- c(4.14, 8.52, 16.31, 32.18, 64.62, 98.76, 151.13, 224.74, 341.35, 
       423.36, 522.78, 674.32, 782.04, 920.01)
DAT <- data.frame(x, y)
mod4 <- onls(y ~ b1 * 10^(b2 * x/(b3 + x)), data = DAT, 
             start = list(b1 = 1, b2 = 5, b3 = 100))


###################################################
### code chunk number 4: onls.Rnw:95-96
###################################################
coef(mod4)


###################################################
### code chunk number 5: onls.Rnw:100-101
###################################################
deviance_o(mod4) 


###################################################
### code chunk number 6: onls.Rnw:105-108
###################################################
mod5 <- onls(y ~ b1 * 10^(b2 * x/(b3 + x)), data = DAT, 
             start = list(b1 = 1, b2 = 5, b3 = 100), extend = c(0, 0))
check_o(mod5)


###################################################
### code chunk number 7: onls.Rnw:111-118
###################################################
x <- 1:100
y <- x^2
set.seed(123)
y <- sapply(y, function(a) rnorm(1, a, 0.1 * a))
DAT <- data.frame(x, y)
mod6 <- onls(y ~ x^a, data = DAT, start = list(a = 1))
mod6


###################################################
### code chunk number 8: onls.Rnw:121-123
###################################################
mod7 <- onls(y ~ x^a, data = DAT, start = list(a = 10), window = 17)
mod7


###################################################
### code chunk number 9: onls.Rnw:128-129
###################################################
plot(mod1)


###################################################
### code chunk number 10: onls.Rnw:132-133
###################################################
plot(mod1, xlim = c(0, 1), ylim = c(0, 1), asp = 1)


###################################################
### code chunk number 11: onls.Rnw:137-138
###################################################
summary(mod1)


###################################################
### code chunk number 12: onls.Rnw:141-142
###################################################
coef(mod1)


###################################################
### code chunk number 13: onls.Rnw:145-146
###################################################
vcov(mod1)


###################################################
### code chunk number 14: onls.Rnw:149-150
###################################################
predict(mod1, newdata = data.frame(conc = 6))


###################################################
### code chunk number 15: onls.Rnw:153-154
###################################################
confint(mod1)


###################################################
### code chunk number 16: onls.Rnw:160-161
###################################################
residuals(mod1)


###################################################
### code chunk number 17: onls.Rnw:164-165
###################################################
fitted(mod1)


###################################################
### code chunk number 18: onls.Rnw:168-169
###################################################
deviance(mod1)


###################################################
### code chunk number 19: onls.Rnw:172-173
###################################################
logLik(mod1)


###################################################
### code chunk number 20: onls.Rnw:179-180
###################################################
residuals_o(mod1)


###################################################
### code chunk number 21: onls.Rnw:183-184
###################################################
deviance_o(mod1)


###################################################
### code chunk number 22: onls.Rnw:187-188
###################################################
logLik_o(mod1)


###################################################
### code chunk number 23: onls.Rnw:194-195
###################################################
x0(mod1)


###################################################
### code chunk number 24: onls.Rnw:198-199
###################################################
y0(mod1)


###################################################
### code chunk number 25: onls.Rnw:205-206
###################################################
plot(mod1)


###################################################
### code chunk number 26: onls.Rnw:212-213
###################################################
plot(mod1, xlim = c(0, 0.25), ylim = c(0, 0.25), asp = 1)


