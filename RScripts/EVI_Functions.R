gauss.curve <- function(parms, tab) {
  t <- 1:nrow(tab)
  parms <- as.list(parms)
  fit1 <- 1 - exp(-((parms$a1 - t[1:(which(t==floor(parms$a1)))])/parms$a4)^parms$a5)
  fit2 <- 1 - exp(-((t[which(t==floor(parms$a1)):length(t)]-parms$a1)/parms$a2)^parms$a3)
  c(fit1, fit2[-1])
}

sigmoid <- function(x, lower_asymptote, carrying_capacity, growth_rate, time_max) {
  return(lower_asymptote + ((carrying_capacity - lower_asymptote)/(1 + exp(-growth_rate * (x - time_max)))))
}


sigmoid.loglik <- function(a, b, c, d) {
  fit <- sigmoid(1:nrow(tabInter), a, b, c, d)
  -sum(dnorm(tabFit$y, rep(fit, ncol(tabInter[,-c(1,2)])), rep(0.25, ncol(tabInter[,-c(1,2)])), log = T), na.rm=T)
} 