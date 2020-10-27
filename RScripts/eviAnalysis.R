library(raster)
library(rgeos)
library(rgdal)
library(geosphere)
library(snow)
library(zoo)
library(circular)
library(biwavelet)
source('RScripts/EVI_Functions.R')

load("Results/eviRaw_4km_2004_2020.RData")
load("Results/snowRaw_4km_2004_2020.RData")

eviRaw$evi[eviRaw$evi<0] <- NA
# matplot(t(eviRaw$evi[1:100,]), type = "o", pch = 16, lwd = 1, lty = 1, col = "darkgreen")

### distance matrix
distM <- geosphere::distm(eviRaw$crds)

res <- do.call("rbind", pbmcapply::pbmclapply(1:nrow(distM), function(i) { 

  ind    <- which(distM[i,]/1000 < 10)
  weight <- approx(c(0,10), c(1,0), distM[i,ind]/1000)$y
  
  x1     <- eviRaw$evi[ind,]; x1[x1<0.1 & !is.na(x1)] <- 0.1
  tmInd  <- sapply(as.numeric(eviRaw$dates), function(x) which.min(abs(x-as.numeric(snowRaw$dates))))
  s      <- snowRaw$snow[ind,tmInd]
  
  # opar <- par(mar = c(6,6,1,6), las = 1)
  # matplot(eviRaw$dates, t(x1), pch = 16, type = "o", col = "darkgreen", xlab = "", ylab = "")
  # par(new=TRUE)
  # matplot(snowRaw$dates[tmInd], t(s), xlim = range(eviRaw$dates), type = "o", pch = 16, col = "cornflowerblue", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  # axis(4)
  
  x <- t(apply(cbind(x1, s), 1, function(y) {
    eviT <- y[1:ncol(x1)]
    snoT <- y[(ncol(x1)+1):length(y)]
    ifelse(snoT==1, NA, ifelse(snoT%in%c(3:4), 0.1, eviT))
  }))
  
  # par(new=TRUE)
  # matplot(eviRaw$dates, t(x), pch = 16, type = "o", col = "orange", xlab = "", ylab = "")
  # par(opar)
  # plot(x[1,], type = "o")

  if(sum(is.na(x[1,]))<(ncol(x)/4) & length(ind)>2 & 
     median(apply(x,1,function(x) var(x, na.rm = T)), na.rm = T)>0.001 &
     findInterval(diff(quantile(x[1,], probs = c(0.025, 0.975), na.rm = T)), c(0.25, 1))) {
  
    x0 <- loess(y ~ t, data = data.frame(t = rep(1:ncol(x), each = nrow(x)), y= c(x)), 
                                         weight = rep(weight, each = ncol(x)), span = 0.02, na.rm = T)
    x1 <- predict(x0, newdata = data.frame(t = 1:ncol(x)))
    x1[x1<0.1] <- 0.1
    # lines(eviRaw$dates, x1, col = "red")
    
    Mx <- x1 - mean(x1, na.rm = T)
    
    fit0  <- optim(fn = leastS.cos, par = c(a = 25, b = 0), f = 52, Mx = Mx, sd = 0.001)
    curve <- fit0$par[1]*cos(pi*((1:length(Mx))/(length(Mx)/((length(Mx)/52)*2))) +
                               (pi+fit0$par[2])) +  mean(fit$time.series[,1], na.rm=T)
    

    if(max(curve)>0.02) {
      
      spl <- which(diff(curve[-length(curve)])<0 & diff(curve[-1])>0 | 
                     diff(curve[-length(curve)])>0 & diff(curve[-1])<0) +1
      
      
      data <- data.frame(t = 1:ncol(x), 
                         date = eviRaw$dates,
                         cos.fit = curve, 
                         y = x1)
      
      data$split1 <- ifelse(1:nrow(data)%in%spl, 1, 0)
      data$split2 <- cut(data$t, breaks = c(0, spl, nrow(data)), labels = F)
      
      data$rise <- as.vector(unlist(apply(cbind(
        as.vector(unlist(lapply(split(data, f=data$split2), function(z) ifelse(coef(lm(z[,"cos.fit"]~z[,"t"]))[2] > 0, TRUE, FALSE)))),
        unlist(lapply(split(data, f=data$split2), function(k) nrow(k)))), 1, function(r) rep(r[1], r[2]))))
      
      splitL <- split(data[data$rise==1,], data$split2[data$rise==1])
      
      
      phen     <- do.call("rbind", lapply(splitL, function(y) {
        
        # y <- splitL[[1]]
        # plot(y$date, y$y)
        
        tab <-merge(data.frame(date = seq(min(y$date), max(y$date), by = "week"), 
                                    t = 1:length(seq(min(y$date), max(y$date), by = "week"))), 
                         y[,c("date", "y")], all.x = T)
          
        if(sum(is.na(tab$y))<nrow(tab)/3) {
          
        st <- tab$date[min(which(tab$y>=(min(tab$y, na.rm = T)+diff(range(tab$y, na.rm = T))/3)))]
          
        # mle <- suppressWarnings(bbmle::mle2(sigmoid.loglik, method="L-BFGS-B",
        #                                     start = list(a = min(tab$y, na.rm = T), 
        #                                                  b = max(tab$y, na.rm = T), 
        #                                                  c = 1, d = 20)))
        # 
        # 
        # plot(tab$t, tab$y)
        # lines(1:nrow(tab), sigmoid(1:nrow(tab), mle@coef[1], mle@coef[2], mle@coef[3], mle@coef[4]))
        # 
        # c(as.numeric(format(median(tab$date, na.rm = T), "%Y")),
            # as.numeric(format(as.Date(approx(1:nrow(tab), tab$date, mle@coef[4])$y), "%j")))
        c(as.numeric(format(median(tab$date, na.rm = T), "%Y")),
          as.numeric(format(as.Date(st), "%j")))
          
        } else c(as.numeric(format(median(tab$date, na.rm = T), "%Y")), NA)
        }))
      
      phenOut  <- merge(data.frame(year = 2004:2020), data.frame(year = phen[,1], start = phen[,2]))
      
    } else phenOut <- cbind(2004:2020, NA)
  } else phenOut <- cbind(2004:2020, NA)
  
  phenOut[,2]
  
}, mc.cores = parallel::detectCores()))

# matplot(2004:2020, t(res), pch = 16, type = "o", col = "grey80")

envOut <- list(crds = eviRaw$crds, years = 2004:2020, evi = res)
save(envOut, file = "eviStart_4km.RData")
