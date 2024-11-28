files <- list.files(here("data", "data"), pattern = "average", full.names = TRUE)

aqua_fun <- function (maxsst, minsst, maxdepth, mindepth) {
  # read in rasters 
  sst <- c(rast(files))
  wc_eez <- st_read(here("data", "data", 'wc_regions_clean.shp'))
  depth <- rast(here("data", "data", "depth.tif"))
  # reproject
  sst_stack_recrs <- project(sst, crs(depth))
  wc_eez <- st_transform(wc_eez, crs(depth))
  
  # find mean sst
  mean_sst <- app(sst_stack, mean)
  # switch to cel.
  mean_sst_c <- mean_sst - 273.15
  # crop depth
  depth_crop <- crop(depth, mean_sst_c)
  # resample depth
  depth_resample <- resample(depth_crop, mean_sst_c, method = "near")
  # suitable sst matrix
  animal_sst <- matrix(c(-Inf, minsst, 0,
                         minsst, maxsst, 1, 
                         maxsst, Inf, 0),
                       ncol = 3, byrow = TRUE)
  # suitable depth matrix
  animal_depth <- matrix(c(-Inf, mindepth, 0,
                         mindepth, maxdepth, 1, 
                         maxdepth, Inf, 0),
                       ncol = 3, byrow = TRUE)
  # reclassify
  reclassified_sst <- classify(mean_sst_c, rcl = animal_sst)
  reclassified_depth <- classify(depth_resample, rcl = animal_depth)
  # multiple together to combine
  reclass_all <- reclassified_sst * reclassified_depth
  # rasterize eez
  wc_eez_rast <- rasterize(wc_eez, reclass_all, "rgn")
  # crop reclassified data to eez
  cropped_animal_locations <- crop(reclass_all, wc_eez_rast)
  # mask animal locations to only eez
  animal_eez <- mask(cropped_animal_locations, wc_eez_rast)
  # find cell area of animal_eez
  area_cell <- cellSize(animal_eez, unit = "km")
  # find the answer :)
  the_answer <- animal_eez * area_cell
  # plot at the end with species name
  zonal_answer <- zonal(the_answer, wc_eez_rast, "sum", na.rm = TRUE) 
  zonal_answer 
}

# oysters:
#   sea surface temperature: 11-30°C
# depth: 0-70 meters below sea level

# red abalone: 
#   sea surface temperature: 8°C - 18°C 
# depth: 0-24 meters below sea level
(maxsst, minsst, maxdepth, mindepth)
aqua_fun(30, 11, 0, -70)

#   transform crs
#   process SST (mean, Cel)
#   crop, resample
#   reclassify
#   sst depth multiplication
#   rasterize eez
#   crop and mask to eez
#   find area of grid cells
#   find total suitable area


