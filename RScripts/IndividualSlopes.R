



indTab <- read.table("/Users/slisovsk/Desktop/WSC_13yrs_3plus.txt", header = T)

head(indTab)
with(indTab, plot(Year, YrDay, pch = 16, col = adjustcolor("cornflowerblue", alpha.f = 0.25)))

  ### all slope
  mod <- lm(YrDay ~ Year, data = indTab)
  
  nsim <- 1000
  bsim <- arm::sim(mod, n.sim = nsim)
  for(i in 1:nsim) abline(coef(bsim)[i,1], coef(bsim)[i,2], col = rgb(0,0,0,0.05))
  
  newdat    <- data.frame(Year = seq(2008, 2020, by = 0.1))  
  newmodmat <- model.matrix(~Year, data = newdat)   
  ftmat     <- matrix(ncol = nsim, nrow = nrow(newdat))
  for(i in 1:nsim) ftmat[,i] <- newmodmat %*% coef(bsim)[i,]   
  with(indTab, plot(Year, YrDay, pch = 16, col = adjustcolor("grey50", alpha.f = 0.25)))
  lines(newdat$Year, apply(ftmat, 1, quantile, prob = 0.5))
  lines(newdat$Year, apply(ftmat, 1, quantile, prob = 0.025), lty = 2)
  lines(newdat$Year, apply(ftmat, 1, quantile, prob = 0.975), lty = 2)

  newy <- matrix(ncol = nsim, nrow = nrow(newdat))  
  for(i in 1:nsim) newy[,i] <- rnorm(nrow(newy), mean = ftmat[,i], sd = bsim@sigma[i])  
  lines(newdat$Year, apply(newy, 1, quantile, prob = 0.025), lty = 3)  
  lines(newdat$Year, apply(newy, 1, quantile, prob = 0.975), lty = 3)  
  
  
  indTab <- read.table("ConklinData/Conklin_IndividDepart.txt", header = T)
  with(indTab, plot(Year, YrDay, pch = 16, col = adjustcolor("cornflowerblue", alpha.f = 0.25)))
  splt <- split(indTab, indTab$BirdID)
  probs <- do.call("cbind", lapply(splt, function(x) {
    if(nrow(x) > 2){
    mod    <- lm(YrDay~Year, data = x)
    newdat <- data.frame(Year = seq(min(x$Year), max(x$Year), by = 0.1))
    lines(newdat$Year, predict(mod, newdata = newdat), col = adjustcolor("orange", alpha.f = 0.5))
    
    bsim <- arm::sim(mod, n.sim = nsim)
    coef(bsim)[,2]
    }
  }))
  quantile(probs, prob = c(0.025, 0.5, 0.975))
  sum(probs<0)/(nsim*ncol(probs))
  
  hist(apply(probs,2,function(x) sum(x<0))/nsim, breaks = seq(0, 1, by = 0.05), main = "", col = adjustcolor("grey50", alpha.f = 0.25), bty = "o", xlab = "Probability of individual slopes being < 0",
       las = 1, cex.lab = 1.2)

  
  slps <- do.call("rbind", (lapply(splt, function(x) {
    if(nrow(x) > 2){
      mod    <- lm(YrDay~Year, data = x)
      bsim <- arm::sim(mod, n.sim = nsim)
      suppressWarnings(cbind(x[1,], n = nrow(x), beta = coef(bsim)[,2]))
    }
  })))
  slps <- subset(slps, abs(beta)<50)
  
  plot(slps$Bill, slps$beta, pch = 16, col = adjustcolor("grey80", alpha.f = 0.15))
  mod <- lmer(beta ~ Bill + Sex + Bill:Sex + (1|BirdID), data = slps, REML = FALSE)
  round(fixef(mod), 3)
  
  bsim <- arm::sim(mod, n.sim = nsim)
  str(bsim)
  round(apply(bsim@fixef, 2, quantile, prob = c(0.025, 0.5, 0.975)), 3)
  
  
  newdat <- expand.grid(Bill = seq(min(slps$Bill), max(slps$Bill), by = 0.25),
                        Sex  = factor(c("m", "f"), levels = levels(slps$Sex)))  
  Xmat   <- model.matrix(~ Bill + Sex + Bill:Sex, data = newdat)   
  ftmat  <- matrix(ncol = nsim, nrow = nrow(newdat))
  for(i in 1:nsim) ftmat[,i] <- Xmat %*% bsim@fixef[i,]
  
  newdat$lower <- apply(ftmat, 1, quantile, prob = 0.025)
  newdat$upper <- apply(ftmat, 1, quantile, prob = 0.975)
  newdat$fit   <- Xmat %*% fixef(mod)
  
  opar <- par(mfrow = c(2,1), las = 1)
  with(slps[slps$Sex=="m",], plot(Bill, beta, pch  =16, col = adjustcolor("cornflowerblue", alpha.f = 0.01), ylim = c(-5, 5), xlim = range(slps$Bill)), ylab = "Ind. Slope in Departure Date")
  with(newdat[newdat$Sex=="m",], lines(Bill, fit, type = "l", col = "cornflowerblue", ylim = range(newdat[,-c(1,2)])))  
  with(newdat[newdat$Sex=="m",], lines(Bill, upper, type = "l", lty = 2, col = "cornflowerblue"))  
  with(newdat[newdat$Sex=="m",], lines(Bill, lower, type = "l", lty = 2, col = "cornflowerblue"))  
  
  with(slps[slps$Sex=="f",], plot(Bill, beta, pch  =16, col = adjustcolor("firebrick", alpha.f = 0.01), ylim = c(-5, 5), xlim = range(slps$Bill)), ylab = "Ind. Slope in Departure Date")
  with(newdat[newdat$Sex=="f",], lines(Bill, fit, type = "l", col = "firebrick", ylim = range(newdat[,-c(1,2)])))  
  with(newdat[newdat$Sex=="f",], lines(Bill, upper, type = "l", lty = 2, col = "firebrick"))  
  with(newdat[newdat$Sex=="f",], lines(Bill, lower, type = "l", lty = 2, col = "firebrick"))  
  par(opar)
  