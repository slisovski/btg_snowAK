#### Stats
library(nlme)
library(lmerTest)
library(lme4)
library(merTools)
library(mgcv)

#### Snow
load("Results/snowMelt_4km.RData")
load("Results/eviStart_4km.RData")

dat <- data.frame(lon = rep(snow$crds[,1], 17), lat = rep(snow$crds[,2], 17), 
                  year = rep(2004:2020, each = nrow(snow$crds)), 
                  site = rep(1:nrow(snow$crds), 17),
                  snow = c(snow$smM),
                  evi  = c(envOut$evi))


plot(dat$year, dat$snow)
plot(dat$year, dat$evi)



### Overall
lm <- lmer(snow ~ year + (1|site), data = dat[dat$year>=2008,], na.action = na.omit)
summary(lm)
confint(lm, parm = "year")

# mean(iav, na.rm = T)
# sd(iav, na.rm = T)



##### north
lm1 <- lmer(snow ~ year + (1|site), data = dat[dat$year>=2008 & dat$lat>=64 & !is.na(dat$snow),])
summary(lm1)
confint(lm1, parm = "year")

##### south
lm2 <- lmer(snow ~ year + (1|site), data = dat[dat$year>=2008 & dat$lat<64 & !is.na(dat$snow),])
summary(lm2)
confint(lm2, parm = "year")


### periods between 2007 and 2014
lm1.2 <- lmer(snow ~ year + (1|site), data = dat[dat$lat>=64 & !is.na(dat$snow) & dat$year%in%c(2008:2014),])
summary(lm1.2)
confint(lm1.2, parm = "year")

lm2.2 <- lmer(snow ~ year + (1|site), data = dat[dat$lat<64 & !is.na(dat$snow) & dat$year%in%c(2008:2014),])
summary(lm2.2)
confint(lm2.2, parm = "year")


newdat  <- data.frame(year = 2008:2020)
mm      <- model.matrix(~year,newdat)

shortdat <- data.frame(year = 2008:2014)
ss       <- model.matrix(~year,shortdat)

opar <- par(mfcol = c(2,2), mar = c(3,3,1,1))
plot(NA, xlim = c(2008, 2020), ylim = c(100, 200), las = 1, ylab = "", xlab = "")
dats <- cbind(2008:2020, t(apply(snow$smM[envOut$crds[,2]>=64,2004:2020>=2008], 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = T))))
polygon(c(dats[,1], rev(dats[,1])), c(dats[,2], rev(dats[,4])), col = adjustcolor("skyblue", alpha.f = 0.75), border = NA)
lines(dats[,1], dats[,3], lwd = 4)

newdat$y1 <- mm%*%fixef(lm1)
with(newdat, lines(year, y1, col = "black", lty = 2, lwd = 1.5))

shortdat$y1 <- ss%*%fixef(lm1.2)
with(shortdat, lines(year, y1, col = "darkorange", lty = 1, lwd = 4.5))


plot(NA, xlim = c(2008, 2020), ylim = c(100, 200), las = 1, ylab = "", xlab = "")
dats <- cbind(2008:2020, t(apply(snow$smM[envOut$crds[,2]<64,2004:2020>=2008], 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = T))))
polygon(c(dats[,1], rev(dats[,1])), c(dats[,2], rev(dats[,4])), col = adjustcolor("skyblue", alpha.f = 0.75), border = NA)
lines(dats[,1], dats[,3], lwd = 4)

newdat$y2 <- mm%*%fixef(lm2)
with(newdat, lines(year, y2, col = "black", lty = 2, lwd = 1.5))

shortdat$y2 <- ss%*%fixef(lm2.2)
with(shortdat, lines(year, y2, col = "darkorange", lty = 1, lwd = 4.5))


####################
#### NDVI ##########
####################

lm <- lmer(evi ~ year + (1|site), data = dat[dat$year>=2008 & dat$year<=2014,], na.action = na.omit)
summary(lm)
confint(lm, parm = "year")

##### north
lm1 <- lmer(evi ~ year + (1|site), data = dat[dat$year>=2008 & dat$year<=2014 & dat$lat>=64 & !is.na(dat$evi),])
summary(lm1)
confint(lm1, parm = "year")

##### south
lm2 <- lmer(evi ~ year + (1|site), data = dat[dat$year>=2008 & dat$lat<64 & !is.na(dat$evi),])
summary(lm2)
confint(lm2, parm = "year")



### periods between 2008 and 2013
lm1.2 <- lmer(evi ~ year + (1|site), data = dat[dat$lat>=64 & !is.na(dat$evi) & dat$year%in%c(2008:2014),])
summary(lm1.2)
confint(lm1.2, parm = "year")

lm2.2 <- lmer(evi ~ year + (1|site), data = dat[dat$lat<64 & !is.na(dat$evi) & dat$year%in%c(2008:2014),])
summary(lm2.2)
confint(lm2.2, parm = "year")


newdata <- data.frame(year = 2008:2020)
mm       <- model.matrix(~year,newdat)

shortdat <- data.frame(year = 2008:2014)
ss       <- model.matrix(~year,shortdat)

plot(NA, xlim = c(2008, 2020), ylim = c(100, 200), las = 1, ylab = "", xlab = "")
dats <- cbind(2008:2020, t(apply(envOut$evi[envOut$crds[,2]>=64,2004:2020>=2008], 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = T))))
polygon(c(dats[!is.na(dats[,2]),1], rev(dats[!is.na(dats[,2]),1])), c(dats[!is.na(dats[,2]),2], rev(dats[!is.na(dats[,2]),4])), col = adjustcolor("seagreen", alpha.f = 0.75), border = NA)
lines(dats[,1], dats[,3], lwd = 4)

newdat$y1 <- mm%*%fixef(lm1)
with(newdat, lines(year, y1, col = "black", lty = 2, lwd = 1.5))

shortdat$y1 <- ss%*%fixef(lm1.2)
with(shortdat, lines(year, y1, col = "darkorange", lty = 1, lwd = 4.5))


plot(NA, xlim = c(2008, 2020), ylim = c(100, 200), las = 1, ylab = "", xlab = "")
dats <- cbind(2008:2020, t(apply(envOut$evi[envOut$crds[,2]<64,2004:2020>=2008], 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = T))))
polygon(c(dats[!is.na(dats[,2]),1], rev(dats[!is.na(dats[,2]),1])), c(dats[!is.na(dats[,2]),2], rev(dats[!is.na(dats[,2]),4])), col = adjustcolor("seagreen", alpha.f = 0.75), border = NA)
lines(dats[,1], dats[,3], lwd = 4)

newdat$y2 <- mm%*%fixef(lm2)
with(newdat, lines(year, y2, col = "black", lty = 2, lwd = 1.5))

shortdat$y2 <- ss%*%fixef(lm2.2)
with(shortdat, lines(year, y2, col = "darkorange", lty = 1, lwd = 4.5))

