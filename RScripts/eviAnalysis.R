library(raster)
library(rgeos)
library(rgdal)
library(geosphere)
library(foreach)
library(snow)
library(zoo)
library(circular)
library(biwavelet)
source('RScripts/EVI_Functions.R')

load("Results/eviRaw_4km.RData")

eviRaw$evi[eviRaw$evi<0] <- NA
matplot(t(eviRaw$evi[1:100,]), type = "o", pch = 16, lwd = 1, lty = 1, col = "darkgreen")

### distance matrix
distM <- geosphere::distm(eviRaw$crds)


res <- do.call("rbind", pbmclapply(1:nrow(distM), function(i) { 

  ind   <- which(distM[i,]/1000 < 7)
  x1     <- eviRaw$evi[ind,]; x1[x1<0 & is.na(x1)] <- NA
  tmInd <- sapply(as.numeric(eviRaw$dates), function(x) which.min(abs(x-as.numeric(snowRaw$dates))))
  s     <- snowRaw$snow[ind,tmInd]
  
  opar <- par(mar = c(6,6,1,6), las = 1)
  matplot(eviRaw$dates, t(x1), pch = 16, type = "o", col = "darkgreen", xlab = "", ylab = "")
  par(new=TRUE)
  matplot(snowRaw$dates[tmInd], t(s), xlim = range(eviRaw$dates), type = "o", pch = 16, col = "cornflowerblue", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  axis(4)
  
  x <- t(apply(cbind(x1, s), 1, function(y) {
    eviT <- y[1:ncol(x1)]
    snoT <- y[(ncol(x1)+1):length(y)]
    ifelse(snoT==1, NA, ifelse(snoT%in%c(3:4), 0, eviT))
  }))
  
  par(new=TRUE)
  matplot(eviRaw$dates, t(x), pch = 16, type = "o", col = "orange", xlab = "", ylab = "")
  par(opar)

  plot(x[1,], type = "o")

  if(sum(is.na(x[1,]))<(ncol(x)/4) & length(ind)>2 & 
     median(apply(x,1,function(x) var(x, na.rm = T)), na.rm = T)>0.001 &
     findInterval(diff(quantile(x[1,], probs = c(0.025, 0.975), na.rm = T)), c(0.25, 1))) {
  
    # matplot(eviRaw$dates, t(x), pch = 16, type = "o", col = "orange", xlab = "", ylab = "")
    x0 <- loess(y ~ t, data = data.frame(t = 1:ncol(x), y= x[1,]), span = 0.02, na.rm = T, weights = ifelse(x[1,]<0.1, 0.1, 1))
    x1 <- predict(x0, newdata = data.frame(t = 1:ncol(x)))
    # lines(eviRaw$dates, x1, col = "red")
    
    
    fit0  <- optim(fn = leastS.cos, par = c(a = 25, b = 0), f = 52, Mx = x1, sd = 0.001)
    curve <- fit0$par[1]*cos(pi*((1:length(x1))/(length(x1)/((length(x1)/52)*2))) +
                               (pi+fit0$par[2])) +  mean(x1, na.rm=T)
  
     # par(new = T)
     # lines(eviRaw$dates, curve)
    
    if(max(curve)>0.02) {
      
      spl <- which(diff(curve[-length(curve)])<0 & diff(curve[-1])>0 | 
                     diff(curve[-length(curve)])>0 & diff(curve[-1])<0) +1
      
      
      data <- cbind(data.frame(t = 1:ncol(x), 
                               date = eviRaw$dates,
                               cos.fit = curve), t(x))
      
      data$split1 <- ifelse(data$t%in%spl, 1, 0) # minima/maxima
      data$split2 <- cut(data$t, breaks = c(0, spl, nrow(data)), labels = F) # segments
      
      ## indicate rise and set periods
      data$rise <- as.vector(unlist(apply(cbind(
        as.vector(unlist(lapply(split(data, f=data$split2), function(z) ifelse(coef(lm(z[,"cos.fit"]~z[,"t"]))[2] > 0, TRUE, FALSE)))),
        unlist(lapply(split(data, f=data$split2), function(k) nrow(k)))), 1, function(r) rep(r[1], r[2]))))
      
      
      splitL0 <- split(data, data$split2)
      splitL  <- splitL0[sapply(splitL0, function(y) unique(y$rise))==1]
      
      
      phenOut0 <- do.call("rbind", lapply(splitL, function(y) {
        
        # y <- splitL[[8]]
        # matplot(1:nrow(y), y[,!is.na(as.numeric(names(y)))], pch = 16, type = "o", col = "orange", xlab = "", ylab = "")
        
        tabInter <-merge(data.frame(date = seq(min(y$date), max(y$date), by = "week"), 
                                    t = 1:length(seq(min(y$date), max(y$date), by = "week"))), 
                         y[,c(2, suppressWarnings(which(!is.na(as.numeric(names(y))))))], all.x = T)
        
        tabFit  <- data.frame(date = rep(tabInter$date, ncol(tabInter)-2), 
                              x = rep(1:nrow(tabInter), ncol(tabInter)-2),
                              y = unlist(c(tabInter[,-c(1:2)])))
        
        if((sum(is.na(tabFit$y))/nrow(tabFit))>0.3) {
          c(as.numeric(format(median(tabInter$date, na.rm = T), "%Y")), NA)
        } else {
          
          
          mle <- suppressWarnings(bbmle::mle2(sigmoid.loglik, method="L-BFGS-B",
                                              start = list(a = min(tabFit$y, na.rm = T), b = max(tabFit$y, na.rm = T), c = 1, d = 20)))
          
          
          # plot(tabFit$x, tabFit$y)
          # lines(1:nrow(tabInter), sigmoid(1:nrow(tabInter), mle@coef[1], mle@coef[2], mle@coef[3], mle@coef[4]))
        
          c(as.numeric(format(median(tabInter$date, na.rm = T), "%Y")), 
            as.numeric(format(as.Date(approx(1:nrow(tabInter), tabInter$date, mle@coef[4])$y), "%j")))
          
        }
          
        
        
      }))
      phenOut  <- merge(data.frame(year = 2004:2020), data.frame(year = phenOut0[,1], start = phenOut0[,2]))
      
    } else phenOut <- cbind(2004:2020, NA)
  } else phenOut <- cbind(2004:2020, NA)
  
  phenOut[,2]
  
}, mc.cores = detectCores()))


res[is.na(as.numeric(res[]))] <- NA 
res[] <- as.numeric(res[])

matplot(2004:2020, t(res), pch = 16, type = "o", col = "grey80")

envOut <- list(crds = eviRaw$crds, years = 2004:2020, snow = snow$smM, evi = res[,c(1,3:10),])
save(envOut, file = "envOur_4km.RData")