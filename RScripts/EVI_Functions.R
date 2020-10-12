gauss.curve <- function(parms, tab) {
  t <- 1:nrow(tab)
  parms <- as.list(parms)
  fit1 <- 1 - exp(-((parms$a1 - t[1:(which(t==floor(parms$a1)))])/parms$a4)^parms$a5)
  fit2 <- 1 - exp(-((t[which(t==floor(parms$a1)):length(t)]-parms$a1)/parms$a2)^parms$a3)
  c(fit1, fit2[-1])
}


fitGauss <- function(tab) {
  
  gauss.loglik <- function(a1, a2, a3, a4, a5) {
    fit <- gauss.curve(parms = list(a1=a1, a2=a2, a3=a3, a4=a4, a5=a5), tab)  
    fit <- ifelse(fit>0.999, 1-(1e-5), ifelse(fit<0.001, 1e-5, fit))
    # cat(paste(c(a1, a2, a3, a4, a5), sep = "  "), "\n")
    -sum(dbinom(x = tab[,2]*100, size = rep(100, length(fit)), prob = fit, log = TRUE), na.rm=T)
  } 
  
  mle <- suppressWarnings(bbmle::mle2(gauss.loglik, method="L-BFGS-B",
                                      start=list(a1 = 200,     a2 = 40,  a3 = 9,  a4 = 40, a5 = 9),
                                      lower=list(a1 = 120,  a2 = 5,   a3 = 0.5,a4 = 5,  a5 = 0.5),
                                      upper=list(a1 = 240,  a2 = Inf,  a3 = Inf, a4 =  Inf, a5 =  Inf)
  ))
  
  coef(mle)
} 


curve_intersect <- function(curve1, curve2, empirical=TRUE, domain=NULL) {
  if (!empirical & missing(domain)) {
    stop("'domain' must be provided with non-empirical curves")
  }
  
  if (!empirical & (length(domain) != 2 | !is.numeric(domain))) {
    stop("'domain' must be a two-value numeric vector, like c(0, 10)")
  }
  
  if (empirical) {
    # Approximate the functional form of both curves
    curve1_f <- approxfun(curve1$x, curve1$y, rule = 2)
    curve2_f <- approxfun(curve2$x, curve2$y, rule = 2)
    
    # Calculate the intersection of curve 1 and curve 2 along the x-axis
    point_x <- uniroot(function(x) curve1_f(x) - curve2_f(x),
                       c(min(curve1$x), max(curve1$x)))$root
    
    # Find where point_x is in curve 2
    point_y <- curve2_f(point_x)
  } else {
    # Calculate the intersection of curve 1 and curve 2 along the x-axis
    # within the given domain
    point_x <- uniroot(function(x) curve1(x) - curve2(x), domain)$root
    
    # Find where point_x is in curve 2
    point_y <- curve2(point_x)
  }
  
  return(list(x = point_x, y = point_y))
} #####


waveletSeas <- function(data) {
  
  spl1 <- cut(1:nrow(data), breaks = c(0, which(data$split1==1 & data$rise==0), nrow(data)), labels = FALSE)
  per  <- as.numeric(format(aggregate(data$date, by = list(spl1), median)$x, "%Y"))
  
  wt  <- wt(cbind(1:nrow(data), na.approx(data[,!is.na(suppressWarnings(as.numeric(names(data))))])))
  # plot(wt)
  power  <- log2(wt$power.corr)
  
  time   <- wt$t 
  period <- wt$period/52
  
  tmp <- foreach(i = unique(spl1), .combine = "rbind", .export = "pck.periods") %dopar% {
    
    if(length(which(spl1==i))>1) {
      tmp01.pow <- apply(power[ ,which(spl1==i)], 1, median, na.rm = T) 
      
      tmp01.sig  <- apply(wt$signif[ ,which(spl1==i)], 1, median, na.rm = T)
      
      ind <- c(1, which((!is.na(tmp01.pow[-length(tmp01.pow)]) &  is.na(tmp01.pow[-1])) |
                          ( is.na(tmp01.pow[-length(tmp01.pow)]) & !is.na(tmp01.pow[-1]))), length(tmp01.pow))
      
      tmp02 <- split(data.frame(period, tmp01.pow, tmp01.sig), f = cut(1:length(tmp01.pow), breaks = unique(ind)))
      
      
      tmp04 <- do.call(rbind, lapply(tmp02, pck.periods))[,1:3]
      
      if(!is.null(tmp04) & nrow(tmp04)>0) {
        
        tmp04 <- tmp04[order(tmp04[,2], decreasing = T),]
        
        tp <- as.numeric(tmp04[1,])
      } else {
        tp <- cbind(NA, NA, NA)
      }  
      
    } else {
      tp <- cbind(NA, NA, NA)
    }
    tp
  }
  
  tmp[c(1,nrow(tmp)),] <- NA
  
  out <- as.data.frame(cbind(per, tmp[,c(1,3,2)]))
  row.names(out) <- NULL
  names(out)     <- c("year", "period", "wt.sign", "wt.power")
  out  
}

predSeas <- function(data, cutoff = 60, info.periods = 4) {
  
  spl1 <- cut(1:nrow(data), breaks = c(0, which(data$split1==1 & data$rise==0), nrow(data)), labels = FALSE)
  per  <- as.numeric(format(aggregate(data$date, by = list(spl1), median)$x, "%Y"))
  
  prds <- as.data.frame(table(spl1))
  
  ind.prds <- suppressWarnings(cbind(min(per):(max(per)-info.periods), ((min(per)-1)+info.periods):(max(per)-1),
                                     prds[-c(1:info.periods),2]))
  
  fc <- foreach:::foreach(loop=1:nrow(ind.prds), .combine = rbind, .packages = "forecast") %dopar% {
    
    if(ind.prds[loop,3]>30) {
      ind <-  which(format(data$date, "%Y")==ind.prds[loop,2]+1)
      ds  <-  data[format(data$date, "%Y")%in%c(ind.prds[loop,1]:ind.prds[loop,2]), names(data)=="1"]
      
      ts  <- ts(ds, frequency = median(ind.prds[loop,3]))
      ets <- stlf(ts, method = "arima", h = length(ind))
      cbind(ind.prds[loop,2]+1, ind, c(ets$mean))
    } else {
      cbind(rep(NA, ind.prds[loop,3]), NA, NA)
    }
  }
  
  out <- data.frame(data[fc[,2], ], as.matrix(fc))
  
  fcy <- foreach:::foreach(loop=unique(out[!is.na(out$V1),"V1"]), .combine = rbind) %dopar% {
    tmp <- out[out[,"V1"]==loop,]
    if(nrow(tmp) < 50) { 
      cbind(prds = loop, pred = NA)
    } else {
      RSS  <- tmp[,names(data)=="1"] - mean(tmp[,names(data)=="1"], na.rm = T)
      SSE  <- sum((tmp[,names(data)=="1"]-tmp[ncol(tmp)])^2, na.rm = T)
      pred <- (sum(RSS^2, na.rm = T)-SSE)/sum(RSS^2, na.rm = T)
      
      cbind(prds  = loop,
            pred  = ifelse(pred<0,0,pred))
    }
  }
  
  
  out  <- cbind(per, fcy[match(per, fcy[,1]),2])
  
  out
}

phenSeas <- function(data) {
  
  spl1 <- cut(1:nrow(data), breaks = c(0, which(data$split1==1 & data$rise==0), nrow(data)), labels = FALSE)
  per  <- as.numeric(format(aggregate(data$date, by = list(spl1), median)$x, "%Y"))
  
  prds <- as.data.frame(table(spl1))
  
  ## Parameterization
  spl1  <- which(data$split1==1 & data$rise==0)	
  tmp1  <- split(data, f = cut(1:nrow(data), breaks = c(0, spl1, nrow(data))))
  z.ind <- suppressWarnings(which(!is.na(as.numeric(names(tmp1[[1]])))))
  
  tmp_loess <-  do.call("rbind", lapply(tmp1, function(z) {
    
    if((nrow(z)<(median(unlist(lapply(tmp1, nrow)))/3)*2) | 
       sum(is.na(apply(z[,z.ind], 1, mean, na.rm = T)))>(nrow(z)/4) |
       median(apply(z[,z.ind],2,var, na.rm = T), na.rm = T)<0.001 |
       var(z$loess, na.rm = T)<0.001) {
      
      rep(NA, 4) 
      
    } else {
      
      if(any(is.na(z$loess))) z$loess <- na.approx(z$loess)
      m1 <- approxfun(x = seq(min(z$loess, na.rm = T), max(z$loess, na.rm = T), length = 100), y = seq(0, 1, length = 100))
      
      ind0 <- which(m1(z$loess)>c(NA,  m1(z$loess)[-nrow(z)]) &  m1(z$loess)>c(m1(z$loess)[-1], NA))
      ind0 <- ind0[which(z[ind0, "loess"]>mean(z[,"loess"], na.rm = T))]
      if(length(ind0)>1) {
        ind0 <- ind0[order(z[ind0, "loess"])][1:2]						
        ind0 <- sort(ind0)
        z[ind0[1]:ind0[2], "loess"] <- NA
        z[(ind0[1]-1):(ind0[2]+1), "loess"] <- na.approx(z[(ind0[1]-1):(ind0[2]+1), "loess"])
      }
      
      diff1 <- c(diff(m1(z$loess)), NA)
      
      sos2 <- z[which.max(diff1), 1]
      eos2 <- z[which.min(diff1), 1]
      
      f.tm2 <- approxfun(x = z[,"t"], y = z[,"yday"])
      c(f.tm2(sos2), f.tm2(eos2), sos2, eos2)
    }
  }))	
  
  circ.year <- approxfun(x = seq(1, 365, length = 360), y = 1:360)
  
  t1 <- circular(tmp_loess[,1], type = "angles", units = "degrees")
  t2 <- attr(median.circular(t1, na.rm = T), which = "medians")[1]
  if(t2>360) t2 <- 360
  year.day  <- approxfun(x = 1:360, y = seq(1, 365, length = 360))
  
  sos_fix <- year.day(t2)
  
  d_loess <- data.frame(x = 1:nrow(tmp_loess), y1 = tmp_loess[,1], y2 = tmp_loess[,2])
  t1 <- circular(d_loess[,3], type = "angles", units = "degrees")
  t2 <- attr(median.circular(t1, na.rm = T), which = "medians")[1]
  if(t2>360) t2 <- 360
  
  eos_fix <- year.day(t2)

  
  ## Split positive curve
  z.ind1 <- suppressWarnings(which(!is.na(as.numeric(names(tmp1[[1]])))))
  
  # fit1 <- foreach(f1=1:length(tmp1), .packages = c("zoo", "bbmle"), .combine = "rbind",
  #                 .export = c("ads.loglik", "ads.curve")) %dopar% {
  
  clusterExport(cl,c("ads.loglik", "ads.curve", "tmp1", "z.ind1", "sos_fix", "eos_fix"))
  fit1 <- do.call("rbind", clusterApply(cl, 1:length(tmp1), function(f1) {
  
                    
                    z <- tmp1[[f1]]
                    
                    if(nrow(z) < 52 | all(z$loess<0.1) | 
                       sum(is.na(apply(z[, z.ind1], 1, mean, na.rm = T)))>(nrow(z)/4) |
                       (quantile(z[,z.ind1], prob = c(0.95), na.rm = T)-quantile(z[, z.ind1], prob = c(0.05), na.rm = T))<0.05) {
                      
                      cbind(rep(NA, nrow(z)), NA, NA)
                      
                    } else {
                      
                      
                      min.t <- quantile(z[, z.ind1], prob = c(0.05), na.rm = T) 
                      max.t <- quantile(z[, z.ind1], prob = c(0.95), na.rm = T)
                      
                      
                      t.temp <- ((1:nrow(z))-mean(1:nrow(z)))
                      mf     <- approxfun(x = z$yday, y = t.temp)
                      
                      sos0 <- sos_fix
                      if(sos0>360) sos0 <- 360
                      
                      eos0 <- eos_fix
                      if(eos0>360) eos0 <- 360
                      
                      
                      mle <- suppressWarnings(mle2(ads.loglik, method = "L-BFGS-B", 
                                                   start=list(c1=quantile(z[,"loess"], prob = 0.1, na.rm = T), 
                                                              c2=diff(quantile(z[,"loess"], prob = c(0, 1), na.rm = T)), 
                                                              w1=0.3, w2=0.3, mu=mf(sos0), v=mf(eos0)),
                                                   upper = list(c1=quantile(z[,"loess"], prob = 0.6, na.rm = T), 
                                                                c2=(max.t-min.t)+0.25, 
                                                                w1 = 2, w2 = 2, mu = mf(sos0)+3, v=mf(eos0)+3),
                                                   lower = list(c1=quantile(z[,"loess"], prob = 0, na.rm = T)-0.01, 
                                                                c2=diff(quantile(z[,"loess"], prob = c(0, 1), na.rm = T)), 
                                                                w1= 0.1, w2 = 0.1, mu = mf(sos0)-3, v=mf(eos0)-3),	
                                                   data = list(y=z[,z.ind1], sd = 0.01)))
                      
                      
                      # matplot(z$date, z[,z.ind1], col = "grey30", lty = 1, pch = 16, cex = 0.5, type = "o")
                      # lines(z$date, z$loess, lwd = 4, col = "firebrick")
                      # lines(z$date, z$y, lwd = 4, col = "darkgreen")
                      #   t0 <- ads.curve(1:nrow(z), parms = list(c1=quantile(z[,"loess"], prob = 0.1, na.rm = T),
                      #              c2=diff(quantile(z[,"loess"], prob = c(0, 1), na.rm = T)),
                      #              w1=0.3, w2=0.3, mu=mf(sos0), v=mf(eos0)))
                      # t1 <- ads.curve(1:nrow(z), parms = list(c1=coef(mle)[1], c2=coef(mle)[2], w1=coef(mle)[3], w2=coef(mle)[4], mu=coef(mle)[5], v=coef(mle)[6]))
                      # lines(z$date, t0, lwd = 3, col = "black", lty = 2)
                      # lines(z$date, t1, lwd = 2, col = "orange")
                      
                      t2 <- approxfun(y = z[,1], x = ((1:nrow(z))-mean(1:nrow(z))))
                      
                      mu <- rep(t2(coef(mle)[5]), nrow(z))
                      v  <- rep(t2(coef(mle)[6]), nrow(z))
                      
                      cbind(ads.curve(seq(1, nrow(z), by = 1), coef(mle)), mu, v)
                      
                    }
                  }))
  
  colnames(fit1) <- c("f1", "f1.mu", "f1.v")
  data <- cbind(data, fit1)
  
  ## Split negative curve
  spl2  <- which(data$split1==1 & data$rise==1)	
  tmp2  <- split(data, f = cut(1:nrow(data), breaks = c(0, spl2, nrow(data))))
  z.ind2 <- suppressWarnings(which(!is.na(as.numeric(names(tmp2[[1]])))))
  
  

  # fit2 <- foreach(f1=1:length(tmp2), .packages = c("zoo", "bbmle"), .combine = "rbind",
  #                   .export = c("ads.loglik", "ads.curve")) %dopar% {        
                    
  clusterExport(cl,c("tmp2", "z.ind2"))
  fit2 <- do.call("rbind", clusterApply(cl, 1:length(tmp2), function(f1) {
     
                    z <- tmp2[[f1]]
                    
                    if((nrow(z) < 52) | 
                       sum(is.na(apply(z[, z.ind2], 1, mean, na.rm = T)))>(nrow(z)/4) |
                       (quantile(z[,z.ind2], prob = c(0.95), na.rm = T)-quantile(z[, z.ind2], prob = c(0.05), na.rm = T))<0.05) {
                      
                      cbind(rep(NA, nrow(z)), NA, NA)
                      
                    }  else {	
                      
                      t.temp <- ((1:nrow(z))-mean(1:nrow(z)))
                      mf <- approxfun(x = z$yday, y = t.temp)
                      
                      sos0 <- sos_fix
                      if(sos0>360) sos0 <- 360
                      
                      eos0 <- eos_fix
                      if(eos0>360) eos0 <- 360
                      
                      mle <- suppressWarnings(mle2(ads.loglik, method = "L-BFGS-B",  
                                                   start = list(c1=quantile(z$loess, prob = 0.9, na.rm = T), 
                                                                c2=-diff(quantile(z$loess, prob = c(0, 1), na.rm = T)), 
                                                                w1=0.3, w2=0.3, mu=mf(eos0), v=mf(sos0)),
                                                   upper = list(c1=quantile(z$loess, prob = 1, na.rm = T)+0.01, 
                                                                c2=-0.01, 
                                                                w1 = 2, w2 = 2, mu= mf(eos0)+3, v=mf(sos0)+3),
                                                   lower = list(c1=quantile(z$loess, prob = 0.6, na.rm = T), 
                                                                c2=-diff(quantile(z$loess, prob = c(0, 1), na.rm = T)), 
                                                                w1 = 0.1, w2 = 0.1, mu = mf(eos0)-3, v=mf(sos0)-3),	
                                                   data = list(y=z[,z.ind2], sd = 0.01)))
                      
                      # matplot(z$date, z[,z.ind2], col = "grey30", lty = 1, pch = 16, cex = 0.5, type = "o")
                      # lines(z$date, z$loess, lwd = 4, col = "firebrick")
                      # lines(z$date, z$y, lwd = 4, col = "darkgreen")
                      #   t0 <- ads.curve(1:nrow(z), parms = list(c1=quantile(z$loess, prob = 0.9, na.rm = T), 
                      #                                           c2=-diff(quantile(z$loess, prob = c(0, 1), na.rm = T)), 
                      #                                           w1=0.3, w2=0.3, mu=mf(eos0), v=mf(sos0)))
                      # t1 <- ads.curve(1:nrow(z), parms = list(c1=coef(mle)[1], c2=coef(mle)[2], w1=coef(mle)[3], w2=coef(mle)[4], mu=coef(mle)[5], v=coef(mle)[6]))
                      # lines(z$date, t0, lwd = 3, col = "black", lty = 2)
                      # lines(z$date, t1, lwd = 2, col = "orange")
                      
                      
                      t2 <- approxfun(y = z[,1], x = ((1:nrow(z))-mean(1:nrow(z))))
                      
                      mu <- rep(t2(coef(mle)[5]), nrow(z))
                      v  <- rep(t2(coef(mle)[6]), nrow(z))
                      
                      cbind(ads.curve(seq(1, nrow(z), by = 1), coef(mle)), mu, v)
                      
                    }
                  }))
  
  colnames(fit2) <- c("f2", "f2.mu", "f2.v")
  data <- cbind(data, fit2)
  
  
  ## smooth transition between split 1 and split 2
  splEnd1 <- which(c(fit1[,1],NA)>c(NA,fit1[,1]) & 
                     c(fit1[,1], NA)>c(fit1[-1,1],NA, NA))
  splEnd2 <- which(c(fit1[,1],NA)<c(NA,fit1[,1]) & 
                     c(fit1[,1], NA)<c(fit1[-1,1],NA, NA))
  splEnd  <- c(0, splEnd1, splEnd2, nrow(data))
  
  tmp3 <- split(data.frame(t = data$t, fit1 = fit1[,1], fit2 = fit2[,1]), 
                f = cut(1:nrow(data), breaks =splEnd[!duplicated(splEnd)], labels = FALSE))
  
  data$fit <- as.numeric(unlist(parLapply(cl, tmp3, function(z) {
    
    # z = tmp3[[1]]
    # plot(z$t, z$fit1, col = "red", type = "o")
    # lines(z$t, z$fit2, col = "blue", type = "o")
    
    if(all(is.na(c(z$fit1, z$fit2)))) rep(NA, nrow(z)) else {
      
      if(sum(is.na(z$fit1))>(nrow(z)/4)) z$fit2 else {
        if(sum(is.na(z$fit2))>(nrow(z)/4)) z$fit1 else { 
          if(z$fit1[1]<z$fit1[nrow(z)]) {
            
            df <- c(z$fit2-z$fit1)
            df[is.na(df)] <- 0
            z$fit1 +df*seq(1, 0, length = nrow(z))
            # lines(z$t, z$fit1 +df*seq(0, 1, length = nrow(z)), col = "orange", type = "l", lwd = 2)
            
          } else {
            
            df <- c(z$fit2-z$fit1)
            df[is.na(df)] <- 0
            z$fit1 +df*seq(0, 1, length = nrow(z))
            # lines(z$t, z$fit1 +df*seq(0, 1, length = nrow(z)), col = "orange", type = "l", lwd = 2)
          }	  					  
          
        }}}
  })))
  
  # plDat <- data[(52*6):(52*15),]
  # with(plDat, plot(date,  y, type = "o"))
  # with(plDat, lines(date, f1,  type = "l", lwd = 3, col = "cornflowerblue"))
  # with(plDat, lines(date, f2,  type = "l", lwd = 3, col = "firebrick"))
  # with(plDat, lines(date, fit, type = "l", lwd = 4, col = "purple"))
  
  
  ## LSP values
  tmp4 <- split(data, f = cut(1:nrow(data), breaks = c(0, spl1, nrow(data))))					
  z.ind <- suppressWarnings(which(!is.na(as.numeric(names(tmp4[[1]])))))
  tmp5 <- matrix(nrow = 6, ncol = length(tmp4))
  
  
  for(phen in 1:length(tmp4)) {
    
    z <- tmp4[[phen]]
    
    # matplot(z[,z.ind], type = "p", cex = 0.5, pch = 16, col = "grey20")
    #   lines(z$fit, col = "orange", lwd = 2)
    
    if(sum(is.na(z$fit)) > nrow(z)/2 |
       (quantile(z[,z.ind], prob = c(0.95), na.rm = T)-quantile(z[, z.ind], prob = c(0.05), na.rm = T))<0.05) {
      
      if(sum(is.na(apply(z[,z.ind], 1, median, na.rm = T)))<(nrow(z)/3) & 
         (quantile(z[,z.ind], prob = c(0.95), na.rm = T)-quantile(z[, z.ind], prob = c(0.05), na.rm = T))>0.05) {
        min.t <- quantile(z[,z.ind], prob = c(0.025), na.rm = T) 
        max.t <- quantile(z[,z.ind], prob = c(0.975), na.rm = T)
      } else {
        min.t <- NA 
        max.t <- NA
      }
      
      tmp5[,phen] <- c(min.t, max.t, rep(NA, 4))
      
    } else { 
      
      f.tm1 <- approxfun(x = 1:nrow(z), y = z[,"t"])
      f.tm2 <- approxfun(x = z[,"t"], y = z[,"yday"])
      
      # amplitude (max, min)
      if(sum(is.na(apply(z[,z.ind], 1, median, na.rm = T)))<(nrow(z)/3)) {
        min.t <- quantile(z[,z.ind], prob = c(0.025), na.rm = T) 
        max.t <- quantile(z[,z.ind], prob = c(0.975), na.rm = T)
      } else {
        min.t <- NA 
        max.t <- NA
      }
      
      tryCatch({
        
        if(all(is.na(z$fit[z$rise==1]))) {
          sos1.t <- NA
          sos2.t <- NA
        } else {
          sos1.t1 <- unique(z$f1.mu[z$rise==1])
          if(length(sos1.t1)>1) sos1.t1 <- sos1.t1[1]
          sos1.t <- f.tm2(sos1.t1)
        }
        
        if(all(is.na(z$fit[z$rise==0]))) {
          eos1.t <- NA
          eos2.t <- NA
        } else {
          eos1.t1 <- unique(z$f1.v[z$rise==0])
          if(length(eos1.t1)>1) eos1.t1 <- eos1.t1[length(eos1.t1)]
          eos1.t <- f.tm2(eos1.t1)
        }
        
        if(!is.na(sos1.t)) {
          ind001 <- c(1:(which.min(abs(z[,1] - sos1.t1))+20))[which.max(z$fit[c(1:(which.min(abs(z[,1] - sos1.t1))+20))])]
          tmp001.1 <- z[1:ind001,]
          f.tmp002 <- approxfun(x = tmp001.1$fit, y = tmp001.1$t)
          sos2.t <-	f.tm2(f.tmp002(min(z$fit[1:ind001], na.rm = T) + 
                                     ((max(z$fit[1:ind001], na.rm = T)-min(z$fit[1:ind001], na.rm = T))*15)/100))
        } else {
          sos2.t <- NA
        }
        
        if(!is.na(eos1.t)) {
          ind001 <- c((which.min(abs(z[,1] - eos1.t1))-20):nrow(z))[which.max(z$fit[(which.min(abs(z[,1] - eos1.t1))-20):nrow(z)])]
          tmp001.1 <- z[ind001:nrow(z),]
          f.tmp002 <- approxfun(x = tmp001.1$fit, y = tmp001.1$t)
          eos2.t <-	f.tm2(f.tmp002(max(z$fit, na.rm = T) - 
                                     ((max(z$fit, na.rm = T)-min(z$fit[ind001:nrow(z)], na.rm = T))*15)/100))		
        } else {
          eos2.t <- NA
        }
        
      },  error = function(e) {
        sos1.t <<- NA
        sos2.t <<- NA
        eos1.t <<- NA
        eos2.t <<- NA})
      
      tmp5[,phen] <- c(min.t, max.t, sos1.t,  sos2.t, eos1.t, eos2.t)
    }	
  }
  
  cbind(per, t(tmp5))   
}







pck.periods <- function(p) {
  if(sum(!is.na(p[,2]))>1) {
    ind01 <- c(ifelse((p[1, 2]>p[2, 2]), 1, NA), which((p[,2]>c(NA, p[-nrow(p),2]) & c(p[-1,2], NA)<p[,2])),
               ifelse(p[nrow(p),2]>p[(nrow(p)-1),2], nrow(p), NA))
    p[ind01[!is.na(ind01)],]
  } else p[which.max(p[,2]),]
}


leastS.cos <- function(params, f, Mx, sd = 0.001) {
  fit  <- params[1]*cos(pi*((1: length(Mx))/(length(Mx)/((length(Mx)/f)*2))) + (pi+params[2]))
  -sum(dnorm(x= Mx[!is.na(Mx)], mean=fit[!is.na(Mx)], sd=sd, log=TRUE))
}


ads.curve <- function(t, parms) {
  t <- (t-mean(t))
  parms <- as.list(parms)
  (parms$c1) + 
    0.5*(parms$c2)*(tanh(parms$w1*(t - parms$mu)) - 
                      tanh(parms$w2*(t - parms$v)))
}


ads.loglik <- function(c1, c2, w1, w2, mu, v) {
  fit <- ads.curve(1:nrow(y), parms = list(c1=c1, c2=c2, w1=w1, w2=w2, mu=mu, v=v))		
  if(ncol(y)>1) {
    -sum(apply(y, 2, function(x) dnorm(x = x, mean = fit, sd = sd, log=TRUE)), na.rm = T)
  } else {
    -sum(dnorm(x = as.vector(y), mean = fit, sd = sd, log=TRUE), na.rm=T)
  }
}

