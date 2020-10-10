library(sf)
library(raster)
library(R.utils)
library(rgdal)

btg <- st_read("data/btg_breedingAK/btg_breedingSiteAK.shp") %>% st_set_crs(4326)
plot(btg$geometry)


fls.gz <- list.files("/bioing/user/slisovsk/4km/", pattern = ".asc.gz", recursive = T,  full.names = T)
dates  <- as.Date(as.POSIXct(unlist(lapply(strsplit(fls.gz, "ims"), function(x) strsplit(x[[2]], "_4km")))[c(TRUE, FALSE)], format = "%Y%j"))

prj <- "+proj=stere +lat_0=90 +lat_ts=60 +lon_0=-80 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6356257 +units=m +no_defs"


  ### Initialize ----
  tab0 <- readLines(fls.gz[600])
  ind <- unlist(suppressWarnings(lapply(tab0, function(x) is.na(as.numeric(gsub(" ", "", x))))))
  tab <- tab0[-which(ind)]
  
  z = unlist(lapply(tab, function(.line) as.numeric(strsplit(.line, '')[[1]])))

  m <- matrix(z, ncol = 6144, nrow = 6144, byrow = T) 
  m   <- m[nrow(m):1,]
  S01 <- raster(m)
  extent(S01) <- c(-12288000, 12288000, -12288000, 12288000)
  proj4string(S01) <- prj
  
  S02 <- S01; S02[] <- 1:length(S02[])
  plot(S01)
  S01_1 <- crop(S01, as(btg %>% st_transform(CRS(prj)) %>% st_geometry(), "Spatial"))
  plot(S01_1)
  plot(btg %>% st_transform(CRS(prj)) %>% st_geometry(), add = T, col = "transparent")
  crdsP <- coordinates(S02)[unlist(raster::extract(S02, as(btg %>% st_transform(CRS(prj)) %>% st_geometry(), "Spatial"))),]
  crds  <- project(crdsP, proj = prj, inv = T) 
  ###### ----


snowM <- matrix(ncol = length(dates), nrow = nrow(crds))
  
for(i in 1:length(dates)) {
    
    cat(sprintf('\rDate %d of %d',
                i, length(dates)))
    
    tab0 <- readLines(fls.gz[i])
    ind <- unlist(suppressWarnings(parallel::mclapply(tab0, function(x) is.na(as.numeric(gsub(" ", "", x))), mc.cores = 15)))
    tab <- tab0[-which(ind)]
    
    z = unlist(parallel::mclapply(tab, function(.line) as.numeric(strsplit(.line, '')[[1]]), mc.cores = 12))
    m <- matrix(z, ncol = 6144, nrow = 6144, byrow = T)
    m   <- m[nrow(m):1,]
    S01 <- raster(m)

    extent(S01) <- c(-12288000, 12288000, -12288000, 12288000)
    proj4string(S01) <- prj

    snowM[,i] <- unlist(raster::extract(S01, crdsP))
    
}
snowRaw <- list(crds = crds, dates = dates, snow = snowM)
save(snowRaw, file = "Results/snowRaw_4km.RData")