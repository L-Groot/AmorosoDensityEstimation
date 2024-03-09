# Questions:

# (1) GG3 should equal Weibull when d = k, but when I set d = k = 2 and Scale = 3
# Weibull looks different from the GG3 and Amorosos

# (2) GG3 should equal Weibull should equal Exponential when d = k = 1. However, when
# I set d = k = 1 and Scale = 2,the Exponential looks different from all the others.


library(devtools)
library(VGAM)
library(AmoRosoDistrib)

# Source dAmoroso() function I wrote
source_url("https://raw.githubusercontent.com/L-Groot/AmorosoDensityEstimation/main/functions.R")

# Set x range
x <- seq(0, 14, by = 0.01)
## Set scale and shape parameters for 3-par GG
Scale <- 4
d <- 2
k <- 2
## Set mu for Amorosos
mu <- 0

# Split plotting window in 5
par(mfrow = c(1,5))


# 3-par Generalized Gamma (3GG)
plot(x, dgengamma.stacy(x, Scale, d = d, k = k), type = "l",
     col = "blue", ylim = 0:1,
     main = "3-par GG",
     sub = "Purple are 5,10,...,95 percentiles", las = 1, ylab = "")
abline(h = 0, col = "blue", lty = 2)
lines(qgengamma.stacy(seq(0.05, 0.95, by = 0.05), Scale, d = d, k = k),
      dgengamma.stacy(qgengamma.stacy(seq(0.05, 0.95, by = 0.05),
                                      Scale, d = d, k = k),
                      Scale, d = d, k = k), col = "purple", lty = 3, type = "h")
lines(x, pgengamma.stacy(x, Scale, d = d, k = k), col = "orange")
abline(h = 0, lty = 2) 


# Amoroso (with density function from AmoRosoDistr package)
# --> should be equal to 3-par GG
plot(x, dgg4(x, a= Scale, l = k, c = d, mu = mu), type = "l",
     col = "blue", ylim = 0:1,
     main = "Amoroso (dgg4())",
     sub = "Purple are 5,10,...,95 percentiles", las = 1, ylab = "")
abline(h = 0, col = "blue", lty = 2)
lines(qgg4(seq(0.05, 0.95, by = 0.05), a= Scale, l = k, c = d),
      dgg4(qgg4(seq(0.05, 0.95, by = 0.05),
                                      a = Scale, l = k, c = d),
                      a = Scale, l = k, c = d), col = "purple", lty = 3, type = "h")
lines(x, pgg4(x, a = Scale, l = k, c = d), col = "orange")
abline(h = 0, lty = 2) 


# Amoroso (with my density function)
# --> should be equal to the two above
plot(x, dAmoroso(x, a= Scale, l = k, c = d, mu = mu), type = "l",
     col = "blue", ylim = 0:1,
     main = "Amoroso (dAmoroso())",
     sub = "Purple are 5,10,...,95 percentiles", las = 1, ylab = "")


# Weibull
# --> should be equal to all of the above if k = d
plot(x, dweibull(x, shape = d, scale = Scale), type = "l",
     col = "blue", ylim = 0:1,
     main = "Weibull (dweibull())",
     sub = "Purple are 5,10,...,95 percentiles", las = 1, ylab = "")


# Exponential
# --> should be equal to all of the above if k = d = 1
# --> YES
plot(x, dexp(x, rate = 1/Scale), type = "l",
     col = "blue", ylim = 0:1,
     main = "Exponential (dexp())",
     sub = "Purple are 5,10,...,95 percentiles", las = 1, ylab = "")



#dgg4(2, a= Scale, l = k, c = d)
#dAmoroso(2, a= Scale, l = k, c = d, mu = 0)
