library(sf)
library(raster)
library(R.utils)

btg <- st_read("data/btg_breedingAK/btg_breedingSiteAK.shp") %>% st_set_crs(4326)
plot(btg$geometry)


filepath <- "/bioing/user/slisovsk/GIS/4km/"
fls      <- list.files(filepath, pattern = "*.zip", recursive = T, full.names = T)

tifs <- data.frame(path = list.files(filepath, pattern = "tif$", recursive = T, full.names = T),
                   year = substr(list.files(filepath, pattern = "tif$", recursive = T, full.names = F), 9, 12),
                   doy  = substr(list.files(filepath, pattern = "tif$", recursive = T, full.names = F), 13, 15))
tifs <- tifs[order(tifs$year, tifs$doy),]
tifs <- tifs[-c(1,2),]


## indices
btgP   <- btg %>% st_transform(CRS(proj4string(raster(tifs[[1]]))))

indTmp <- raster(tif.fls[[1]]); indTmp[] <- 1:length(indTmp[])
crds   <- coordinates(indTmp)[unlist(raster::extract(indTmp, as(btgP$geometry, "Spatial"))),]

# snowOut <- matrix(nrow = nrow(crds), ncol = nrow(tifs))
# 
# for(i in 1:nrow(tifs)){
#   
#   cat(sprintf("\r%d of %d", i, nrow(tifs)))
#   
#   snowOut[,i] <- unlist(raster::extract(raster(as.character(tifs$path[i])), as(btgP$geometry, "Spatial")))
#   
# }
# 
# snowList <- list(year = tifs$year, doi = tifs$doy, crds = crds, dat = snowOut)
# save(snowList, file = "tmp/snowList.RData")
load("tmp/snowList.RData")

