# Load required libraries
library(ggplot2)

# Define parameters for Amoroso so that Amoroso = Weibull
n <- 1000
a <- 3     # Scale parameter
lambda <- 1  # Shape parameter for Weibull (when c>0)
c <- 2     # Shape parameter for Amoroso
mu <- 0.5    # Location parameter (shape par of Weibull must be >0)

# Define density function for Amoroso
dAmoroso <- function(x, a, lambda, c, mu) {
  c1 <- 1/(gamma(lambda))
  c2 <- abs(c/a)
  c3 <- ((x-mu)/a)^(lambda*c-1)
  c4 <- exp(-((x-mu)/a)^c)
  return(c1*c2*c3*c4)
}

# Define inverse cumulative distribution function (CDF) for Amoroso
pAmoroso <- function(q, a, lambda, c, mu) {
  p <- integrate(dAmoroso, lower = mu, upper = q, a = a, lambda = lambda, c = c, mu = mu)$value
  return(p)
}

# Define function to sample from Amoroso
rAmoroso <- function(n, a, lambda, c, mu) {
  samples <- numeric(n)
  for (i in 1:n) {
    u <- runif(1)
    samples[i] <- uniroot(function(x) pAmoroso(x, a, lambda, c, mu) - u, interval = c(mu, 1000))$root
    #samples[i] <- uniroot(function(x) pAmoroso(x, a, lambda, c, mu) - u, interval = c(mu, mu + 10 * a))$root
  }
  return(samples)
}

# Simulate data from Amoroso distribution using inverse transform sampling
amoroso_samples <- rAmoroso(n, a, lambda, c, mu)

# Simulate data from Weibull with eta(shape) = lambda and beta(scale) = c
#a = 2
#mu = 0.5 # because shape par of Weibull must be >0
weibull_samples <- rweibull(n, shape = mu, scale = a)
max(weibull_samples)

# Divide plot window into 2
par(mfrow=c(2,1))

# Determine the range of values for the x-axis
x_min <- min(c(weibull_samples, amoroso_samples))
x_max <- max(c(weibull_samples, amoroso_samples))

# Plot histogram of simulated data for Weibull distribution
hist(weibull_samples, freq = FALSE, col = "lightblue", main = "Simulated Weibull Data", xlab = "Value", breaks = 20, xlim = c(x_min, x_max))

# Plot histogram of simulated data for Amoroso distribution
hist(amoroso_samples, freq = FALSE, col = "red", main = "Simulated Amoroso Data", xlab = "Value", breaks = 20, xlim = c(x_min, x_max))


#--------------------------------------------------------------------------------
# Plot histogram of simulated data
hist(weibull_samples, freq = FALSE, col = "lightblue", main = "Simulated Weibull Data", xlab = "Value", breaks = 20)

# Plot histogram of simulated data
hist(amoroso_samples, freq = FALSE, col = "red", main = "Simulated Amoroso Data", xlab = "Value", breaks = 20)

# Add density curve
#curve(dAmoroso(x, a, lambda, c, mu), add = TRUE, col = "blue", lwd = 2)
