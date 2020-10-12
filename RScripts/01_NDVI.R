library(raster)
library(rgeos)
library(gdalUtils)
library(rgdal)
library(geosphere)
library(foreach)
library(snow)
library(zoo)
library(circular)
library(biwavelet)
library(RNetCDF)
library(maptools); data(wrld_simpl)

### Breeding range polygons
bs0 <- readShapePoly("F:/Archive/BirdBreedingPolygons/Waders/Wader_2.shp")
ind <- which(bs0$SCINAME=="Limosa lapponica" & bs0$SEASONAL==2)

btg  <- bs0[1,]
proj4string(btg) <- proj4string(wrld_simpl)


plot(btg, col = "firebrick")
plot(wrld0, add = T)
box()


load("snow_4km.RData")
load("snowRaw_4km.RData")

dwd <- "F:/GeoDat/MODIS/VHP_SM_SMN"
fls <- list.files(dwd, pattern = "SMN.tif$")

year <- as.numeric(substring(fls, 17, 20))
week <- as.numeric(substring(fls, 21, 23))

date0 <- cbind(year, week)
date <- as.Date(as.POSIXct(apply(date0, 1, function(x) {
  tm <- seq(as.POSIXct(paste0(x[1], "-01-01")), as.POSIXct(paste0(x[1], "-12-31")), by = "day")
  w  <- which(x[2]==as.numeric(format(tm, "%U")))
  mean(tm[w])
}), origin = "1970-01-01"))

files <- data.frame(Path = list.files(dwd, pattern = "SMN.tif$", full.names = T), Year = year, Week = week, Date = date)

dates <- date[year%in%c(2004:2018)]

# eviM  <- matrix(nrow = nrow(snow$crds), ncol = length(dates))
# 
# for(i in 1:length(dates)) {
# 
#   cat(sprintf('\rRaster %d of %d',
#               i, length(dates)))
#   
#   r0 <- raster(as.character(files[which(date==dates[i]),1]))
#   eviM[,i] <- extract(r0, snow$crds)
# }
# 
# 
# eviRaw <- list(crds = snow$crds, dates = dates, evi = eviM)
# save(eviRaw, file = "eviRaw_4km.RData") 
load("eviRaw_4km.RData")

eviRaw$evi[eviRaw$evi<0] <- NA
matplot(t(eviRaw$evi[1:100,]), type = "o", pch = 16, lwd = 1, lty = 1, col = "darkgreen")

res     <- array(dim = c(nrow(eviRaw$crds), 10, 15))

plot(btg, col = NA)
box()

### makeCluster
cl <- makeCluster(rep("localhost", 16), type = "SOCK")
null <- clusterEvalQ(cl, library(zoo))
null <- clusterEvalQ(cl, library(bbmle))
#####

for(i in 1:nrow(eviRaw$crds)) {
 
   
  percent <- i/(nrow(eviRaw$crds))*100
  cat(sprintf('\r[%-50s] %d%% - %d',
              paste(rep('=', percent/2), collapse = ''),
              floor(percent), i))
  
  if(all(is.na(res[i,7,]))) {
  
  ind0 <- apply(eviRaw$crds, 1, function(x) distVincentySphere(x, eviRaw$crds[i,])/1000)
  ind  <- order(ind0)[ind0[order(ind0)]<10]
  
  x <- eviRaw$evi[ind,]; x[x<0 & is.na(x)] <- NA
    tmInd <- apply(cbind(as.POSIXct(eviRaw$dates, tz = "GMT")), 1, function(x) which.min(abs(x-as.numeric(as.POSIXct(snowRaw$dates, tz = "GMT")))))
  s <- snowRaw$snow[ind,tmInd] ## 1 = SEA ## 2 = Land without snow ## 3 = Sea Ice ## 4 = Snow covered land
  
  # opar <- par(mar = c(6,6,1,6), las = 1)
  # matplot(eviRaw$dates, t(x), pch = 16, type = "o", col = "darkgreen", xlab = "", ylab = "")
  # par(new=TRUE)
  # matplot(snowRaw$dates[tmInd], t(s), xlim = range(eviRaw$dates), type = "o", pch = 16, col = "cornflowerblue", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  # axis(4)

  
  x <- t(apply(cbind(x, s), 1, function(y) {
    eviT <- y[1:ncol(x)]
    snoT <- y[(ncol(x)+1):length(y)]
    ifelse(snoT==1, NA, ifelse(snoT%in%c(3:4), 0, eviT))
  }))
  
  # par(new=TRUE)
  # matplot(eviRaw$dates, t(x), pch = 16, type = "o", col = "orange", xlab = "", ylab = "")
  # par(opar)
  # plot(x[1,], type = "o")
  
  if(sum(is.na(x[1,]))<(ncol(x)/4) & length(ind)>2 & 
     median(apply(x,1,function(x) var(x, na.rm = T)), na.rm = T)>0.001 &
     findInterval(diff(quantile(x[1,], probs = c(0.025, 0.975), na.rm = T)), c(0.25, 1))) {
        
        
        x0 <- loess(y ~ t, data = data.frame(t = 1:ncol(x), y= x[1,]), span = 0.02, na.rm = T)
        x1 <- predict(x0, newdata = data.frame(t = 1:ncol(x)))
        # lines(eviRaw$dates, x1, col = "red")
        
        tmp <- data.frame(date = eviRaw$dates, y = na.approx(x[1,]))
        
        ts <- ts(tmp$y, frequency = 52,
                 start = as.numeric(unlist(strsplit(format(tmp$date[1], "%Y %j"), " "))))
        fit = stl(ts, s.window='periodic')
        
        
        Mx <- fit$time.series[,1] - mean(fit$time.series[,1], na.rm = T)
        
        fit0  <- optim(fn = leastS.cos, par = c(a = 25, b = 0), f = 52, Mx = Mx, sd = 0.001)
        curve <- fit0$par[1]*cos(pi*((1:length(Mx))/(length(Mx)/((length(Mx)/52)*2))) +
                                   (pi+fit0$par[2])) +  mean(fit$time.series[,1], na.rm=T)
        
        # par(new=TRUE)
        # plot(eviRaw$dates, curve, col = "red", type = "l")
        
        if(max(curve)>0.02) {
          
          spl <- which(diff(curve[-length(curve)])<0 & diff(curve[-1])>0 | 
                         diff(curve[-length(curve)])>0 & diff(curve[-1])<0) +1
          
          data <- cbind(data.frame(t = 1:ncol(x), date = eviRaw$dates, yday = as.numeric(format(eviRaw$dates, "%j")), loess = x1, cos.fit = curve), t(x))
          
          data$split1 <- ifelse(1:nrow(data)%in%spl, 1, 0) # minima/maxima
          data$split2 <- cut(data$t, breaks = c(0, spl, nrow(data)), labels = F) # segments
          
          ## indicate rise and set periods
          data$rise <- as.vector(unlist(apply(cbind(
            as.vector(unlist(lapply(split(data, f=data$split2), function(z) ifelse(coef(lm(z[,"cos.fit"]~z[,"t"]))[2] > 0, TRUE, FALSE)))),
            unlist(lapply(split(data, f=data$split2), function(k) nrow(k)))), 1, function(r) rep(r[1], r[2]))))
          
          
          data <- cbind(data,
                        data.frame(
                          seasonal  = fit$time.series[,1],
                          trend     = fit$time.series[,2],
                          remainder = fit$time.series[,3]))
          
          seas   <- waveletSeas(data)
          seas   <- cbind(data.frame(year = 2004:2018), seas[match(2004:2018, seas[,1]),-1])
          seas2  <- predSeas(data)
          seas   <- cbind(seas, pred = seas2[match(seas[,1], seas2[,1]),2])  
          
          seas3  <- tryCatch(phenSeas(data), error = function(x) matrix(rep(NA, 7), ncol = 7))
          
          colnames(seas3) <- c("year", "min", "max", "sos.1", "sos.2", "eos.1", "eos.2")
          seas <- cbind(seas, seas3[match(seas[,1], seas3[,1]),-1])  
          
          
          res[i,,] <- t(seas[,-1]) ## [cells, c("1:period", "2:wt.sign", "3:wt.power", "4:pred", "5:min", "6:max", "7:sos.1", "8:sos.2", "9:eos.1", "10:eos.2")]
        }
        points(eviRaw$crds[i,1], eviRaw$crds[i,2], pch = 16, cex = 0.2, col = "firebrick")
    } else points(eviRaw$crds[i,1], eviRaw$crds[i,2], pch = 16, cex = 0.2, col = "grey10")
  
  }

  
  
  if(floor(i/500)==i/500) {
    r0 <- raster(extent(btg), res = 0.1)
      proj4string(r0) <- proj4string(wrld_simpl)
    test <- rasterize(eviRaw$crds, r0, field = res[,1,12], fun = function(x, na.rm) median(x, na.rm = TRUE))
    plot(test)
    plot(btg, add = T)
    save(res, file = "eviAK_4km.RData")
  }
  
}

stopCluster(cl)


envOut <- list(crds = eviRaw$crds, years = 2004:2018, snow = snow$smM, evi = res[,c(1,3:10),])
save(envOut, file = "envOur_4km.RData")










