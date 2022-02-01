load("mu_APOEnb4.RData")
load("sigma_APOEnb4.RData")

time_breaks <- c(40, 65, 70, 75, 80, 85, 90, 95)
time_breaks0 <- c(0, time_breaks)
time_breaksInf <- c(time_breaks, Inf)
time_intervals <- time_breaks0[-1] - time_breaks0[-9]


CumulLambda.APOE <- function(nb4, t, u) {
  if (nb4 == 0) {
    alpha <- exp(mu[,1] + u*sigma[,1])
  } else if (nb4 == 1) {
    alpha <-  exp(mu[,2] + u*sigma[,2])
  } else if (nb4 == 2) {
    alpha <-  exp(mu[,3] + u*sigma[,3])
  } 
  sum1 <- sum((alpha[-9]*time_intervals)[t >= time_breaks])
  sum2 <- ((alpha*(t - time_breaks0))[t < time_breaksInf])[1]
  return(sum1 + sum2)
}

lambda.APOE <- function(nb4, t, u) {
  if (nb4 == 0) {
    alpha <- exp(mu[,1] + u*sigma[,1])
  } else if (nb4 == 1) {
    alpha <-  exp(mu[,2] + u*sigma[,2])
  } else if (nb4 == 2) {
    alpha <-  exp(mu[,3] + u*sigma[,3])
  } 
  return(alpha[max(which(t - time_breaks0 >= 0))])
}

