load("tmp/snowList.RData")

######## Functions -----
library(bbmle)
library(zoo)
library(parallel)
library(doParallel)

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
} 


######## Curve fitting and phenology ------
smM       <- matrix(nrow = nrow(snowList$crds), ncol = length(unique(snowList$year)))

for(j in 1:nrow(smM)) {
  
  cat(sprintf('\rLocation %d of %d',
              j, nrow(smM)))
  
  y <- ifelse(snowList$dat[j,]%in%c(4,165), 4, ifelse(snowList$dat[j,]%in%c(3,164), 3, snowList$dat[j,]))
  
  if(sum(!y%in%c(2,4)) < length(y)/5) {
    
    tab0 <- data.frame(year = snowList$year, 
                       s    = ifelse(y%in%c(0,1), 0, ifelse(y==2, 0, 1)), 
                       doy = snowList$doi)
    
    spl <- split(tab0, f = as.character(tab0$year))
    
    sm1 <- do.call("rbind", mclapply(spl, function(x) {
      
      tab <- merge(data.frame(day = 1:365), data.frame(day = as.numeric(x[,3]), p = as.numeric(x[,2])), all.x = T)
      
      mle <- fitGauss(tab)
      fit <- gauss.curve(mle, tab)
      
      sm <- tryCatch(curve_intersect(data.frame(x = tab[,1], y = fit)[1:mle[1],], data.frame(x = tab[,1], y = 0.666))$x,
                     error = function(x) NA)
      
      cbind(year = median(as.numeric(as.numeric(as.character(x[,1])))), sm = sm)
      
    }, mc.cores = 15))
    
    smM[j,] <- merge(data.frame(year = 2006:2020), as.data.frame(sm1), all.x = T)$sm
    
  }
  
}

snow <- list(crds = crds, smM = smM)
save(snow, file = "/home/slisovsk/Documents/btg_snowAK/tmp/snow_4km.RData")
