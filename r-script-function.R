# List files to read in
files <- list.files(here("data"), pattern = "average", full.names = TRUE)

# Make a function
aqua_fun <- function (maxsst, minsst, maxdepth, mindepth, species) {
  # Read in rasters 
  sst <- c(rast(files))
  wc_eez <- st_read(here("data", 'wc_regions_clean.shp'))
  depth <- rast(here("data", "depth.tif"))
  #Find mean sst
  mean_sst <- app(sst_stack, mean)
  # Switch to cel.
  mean_sst_c <- mean_sst - 273.15
  # crop depth
  depth_crop <- crop(depth, mean_sst_c)
  # Resample depth
  depth_resample <- resample(depth_crop, mean_sst_c, method = "near")
  # Suitable sst matrix
  animal_sst <- matrix(c(-Inf, minsst, 0,
                         minsst, maxsst, 1, 
                         maxsst, Inf, 0),
                       ncol = 3, byrow = TRUE)
  # Suitable depth matrix
  animal_depth <- matrix(c(-Inf, mindepth, 0,
                           mindepth, maxdepth, 1, 
                           maxdepth, Inf, 0),
                         ncol = 3, byrow = TRUE)
  # Reclassify
  reclassified_sst <- classify(mean_sst_c, rcl = animal_sst)
  reclassified_depth <- classify(depth_resample, rcl = animal_depth)
  reclassified_sst[reclassified_sst == 0] <- NA
  reclassified_depth[reclassified_depth == 0] <- NA
  # rematch CRSs
  reclassified_depth <- project(reclassified_depth, reclassified_sst)
  # Multiple together to combine
  reclass_all <- reclassified_sst * reclassified_depth
  # Rasterize eez
  wc_eez_rast <- rasterize(wc_eez, reclass_all, "rgn")
  # Crop reclassified data to eez
  cropped_animal_locations <- crop(reclass_all, wc_eez_rast)
  # Mask animal locations to only eez
  animal_eez <- mask(cropped_animal_locations, wc_eez_rast)
  # Find cell area of animal_eez
  area_cell <- cellSize(animal_eez, 
                        mask = TRUE,
                        unit = "km")
  # Multiply oyster data by area_cell data to find area of eez
  eez_area <- animal_eez * area_cell
  
  # Use zonal function to separate area by EEZ
  zonal_area <- zonal(eez_area, wc_eez_rast, "sum", na.rm = TRUE) 
  
  # Print table
  kable(eez_suitable, digits = 2,
        caption = "Suitable area in West Coast EEZs", 
        col.names = c("EEZ Region", "Area (km^2)"))
  
  # Join for plotting
  eez <- left_join(wc_eez, eez_suitable, by = "rgn")
  
  # Map suitable locations by EEZ
  tm_shape(depth) +
    tm_raster(palette = "-GnBu",
              title = "Bathymetry\n(m above and below sea level)",
              midpoint = 0,
              legend.show = FALSE) +
    tm_shape(eez, raster.downsample = TRUE) +
    tm_polygons(col = "mean",
                palette = "Purples",
                alpha = 0.8,
                linewidth = 0.2,
                title = expression("Suitable habitat area (km"^2*")")) +
    tm_text("rgn", size = 0.45) +
    tm_compass(size = .5,
               position = c("left", "bottom")) +
    tm_scale_bar(position = c("left", "bottom")) +
    tm_layout(legend.outside = TRUE,
              frame = FALSE,
              main.title = paste0("Suitable Habitats by West Coast EEZs\n(Exclusive Economic Zones) for ", species),
              main.title.size = .7)
}
