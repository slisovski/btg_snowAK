library(sf)

btg <- st_read("data/btg_breedingAK/btg_breedingSiteAK.shp") %>% st_set_crs(4326)
plot(btg$geometry)


fielpath <- ""


