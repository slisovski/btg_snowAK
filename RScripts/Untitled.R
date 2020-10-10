#### Stats
library(nlme)
library(lmerTest)
library(lme4)
library(merTools)
library(mgcv)

#### Snow
load("Results/snowMelt_4km.RData")

dat <- data.frame(lon = rep(snow$crds[,1], 13), lat = rep(snow$crds[,2], 13), 
                  year = rep(2008:2020, each = nrow(snow$smM[,-c(1,2)])), 
                  site = rep(1:nrow(snow$crds), 13),
                  snow = c(snow$smM[,-c(1,2)]))


plot(dat$year, dat$snow)


### Overall
lm <- lmer(snow ~ year + (1|site), data = dat, na.action = na.omit)
summary(lm)
confint(lm, parm = "year")

# mean(iav, na.rm = T)
# sd(iav, na.rm = T)



##### north
lm1 <- lmer(snow ~ year + (1|site), data = dat[dat$lat>=64 & !is.na(dat$snow),])
summary(lm1)
confint(lm1, parm = "year")

##### south
lm2 <- lmer(snow ~ year + (1|site), data = dat[dat$lat<64 & !is.na(dat$snow),])
summary(lm2)
confint(lm2, parm = "year")


### periods between 2007 and 2013
lm1.2 <- lmer(snow ~ year + (1|site), data = dat[dat$lat>=64 & !is.na(dat$snow) & dat$year%in%c(2007:2013),])
summary(lm1.2)
confint(lm1.2, parm = "year")

lm2.2 <- lmer(snow ~ year + (1|site), data = dat[dat$lat<64 & !is.na(dat$snow) & dat$year%in%c(2007:2013),])
summary(lm2.2)
confint(lm2.2, parm = "year")


newdat  <- data.frame(year = 2004:2018)
mm      <- model.matrix(~year,newdat)

shortdat <- data.frame(year = 2007:2013)
ss       <- model.matrix(~year,shortdat)

opar <- par(mfcol = c(2,2), mar = c(3,3,1,1))
plot(NA, xlim = c(2004, 2018), ylim = c(100, 200), las = 1, ylab = "", xlab = "")
dats <- cbind(2004:2018, t(apply(envOut$snow[envOut$crds[,2]>=64,], 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = T))))
polygon(c(dats[,1], rev(dats[,1])), c(dats[,2], rev(dats[,4])), col = adjustcolor("skyblue", alpha.f = 0.75), border = NA)
lines(dats[,1], dats[,3], lwd = 4)

newdat$y1 <- mm%*%fixef(lm1)
with(newdat, lines(year, y1, col = "black", lty = 2, lwd = 1.5))

shortdat$y1 <- ss%*%fixef(lm1.2)
with(shortdat, lines(year, y1, col = "darkorange", lty = 1, lwd = 4.5))


plot(NA, xlim = c(2004, 2018), ylim = c(100, 200), las = 1, ylab = "", xlab = "")
dats <- cbind(2004:2018, t(apply(envOut$snow[envOut$crds[,2]<64,], 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = T))))
polygon(c(dats[,1], rev(dats[,1])), c(dats[,2], rev(dats[,4])), col = adjustcolor("skyblue", alpha.f = 0.75), border = NA)
lines(dats[,1], dats[,3], lwd = 4)

newdat$y2 <- mm%*%fixef(lm2)
with(newdat, lines(year, y2, col = "black", lty = 2, lwd = 1.5))

shortdat$y2 <- ss%*%fixef(lm2.2)
with(shortdat, lines(year, y2, col = "darkorange", lty = 1, lwd = 4.5))


####################
#### NDVI ##########
####################


##### north
lm1 <- lmer(eviSos1 ~ year + (1|site), data = dat[dat$lat>=64 & !is.na(dat$eviSos1),])
summary(lm1)
confint(lm1, parm = "year")

##### south
lm2 <- lmer(eviSos1 ~ year + (1|site), data = dat[dat$lat<64 & !is.na(dat$eviSos1),])
summary(lm2)
confint(lm2, parm = "year")


### periods between 2008 and 2013
lm1.2 <- lmer(eviSos1 ~ year + (1|site), data = dat[dat$lat>=64 & !is.na(dat$eviSos1) & dat$year%in%c(2007:2013),])
summary(lm1.2)
confint(lm1.2, parm = "year")

lm2.2 <- lmer(eviSos1 ~ year + (1|site), data = dat[dat$lat<64 & !is.na(dat$eviSos1) & dat$year%in%c(2007:2013),])
summary(lm2.2)
confint(lm2.2, parm = "year")


newdata <- data.frame(year = 2004:2018)
mm       <- model.matrix(~year,newdat)

shortdat <- data.frame(year = 2007:2013)
ss       <- model.matrix(~year,shortdat)

plot(NA, xlim = c(2004, 2018), ylim = c(100, 200), las = 1, ylab = "", xlab = "")
dats <- cbind(2004:2018, t(apply(envOut$evi[envOut$crds[,2]>=64,6,], 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = T))))
polygon(c(dats[!is.na(dats[,2]),1], rev(dats[!is.na(dats[,2]),1])), c(dats[!is.na(dats[,2]),2], rev(dats[!is.na(dats[,2]),4])), col = adjustcolor("seagreen", alpha.f = 0.75), border = NA)
lines(dats[,1], dats[,3], lwd = 4)

newdat$y1 <- mm%*%fixef(lm1)
with(newdat, lines(year, y1, col = "black", lty = 2, lwd = 1.5))

shortdat$y1 <- ss%*%fixef(lm1.2)
with(shortdat, lines(year, y1, col = "darkorange", lty = 1, lwd = 4.5))


plot(NA, xlim = c(2004, 2018), ylim = c(100, 200), las = 1, ylab = "", xlab = "")
dats <- cbind(2004:2018, t(apply(envOut$evi[envOut$crds[,2]<64,6,], 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = T))))
polygon(c(dats[!is.na(dats[,2]),1], rev(dats[!is.na(dats[,2]),1])), c(dats[!is.na(dats[,2]),2], rev(dats[!is.na(dats[,2]),4])), col = adjustcolor("seagreen", alpha.f = 0.75), border = NA)
lines(dats[,1], dats[,3], lwd = 4)

newdat$y2 <- mm%*%fixef(lm2)
with(newdat, lines(year, y2, col = "black", lty = 2, lwd = 1.5))

shortdat$y2 <- ss%*%fixef(lm2.2)
with(shortdat, lines(year, y2, col = "darkorange", lty = 1, lwd = 4.5))




lm <- lmer(snow ~ eviSos1 + (1|site), data = dat, na.action = na.omit)
summary(lm)
MuMIn::r.squaredGLMM(lm)


res <- cor.test(dat$snow, dat$eviSos1, method = "pearson")
res

resN <- cor.test(dat$snow[dat$lat>=64], dat$eviSos1[dat$lat>=64], method = "pearson")
resN

resS <- cor.test(dat$snow[dat$lat<64], dat$eviSos1[dat$lat<64], method = "pearson")
resS
