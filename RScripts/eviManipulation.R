library(raster)
library(rgeos)
library(gdalUtils)
library(rgdal)
library(sf)
library(parallel)
library(RNetCDF)
library(pbmcapply)  

### Breeding range polygons
btg <- st_read("data/btg_breedingAK/btg_breedingSiteAK.shp") %>% st_set_crs(4326)

load("Results/snowRaw_4km_2004_2020.RData")

dwd <- "/bioing/user/slisovsk/VHP_SM_SMN/"
fls <- list.files(dwd, pattern = "SMN.tif$")

strs <- sapply(strsplit(fls, ".", fixed = T), function(x) x[[5]])


year <- as.numeric(substring(strs, 2, 5))
week <- as.numeric(substring(strs, 6, 8))

date0 <- cbind(year, week)
date <- as.Date(as.POSIXct(apply(date0, 1, function(x) {
  tm <- seq(as.POSIXct(paste0(x[1], "-01-01")), as.POSIXct(paste0(x[1], "-12-31")), by = "day")
  w  <- which(x[2]==as.numeric(format(tm, "%U")))
  mean(tm[w])
}), origin = "1970-01-01"))

files <- data.frame(Path = list.files(dwd, pattern = "SMN.tif$", full.names = T), Year = year, Week = week, Date = date)
files <- subset(files, Year >= 2004)
files <- files[order(files$Date),]


eviM <- do.call("cbind", pbmclapply(1:nrow(files), function(x) { 
   r0 <- raster(as.character(files[x,1]))
   raster::extract(r0, snowRaw$crds)
}, mc.cores = detectCores()-10))

eviRaw <- list(crds = snow$crds, dates = dates, evi = eviM)
save(eviRaw, file = "Results/eviRaw_4km.RData") 