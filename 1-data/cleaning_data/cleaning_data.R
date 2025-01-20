library(dplyr)
library(ggplot2)
library(sf)
library(terra)
library(raster)
library(mapview)
library(tidyr)
library(rstudioapi)
# Set working directory
data_path <- normalizePath(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/data/"))
setwd(data_path)

# Generate a raster data
study_area <- read_sf('Boundary_study_area.shp') # study area

crs_wgs84 <- CRS("+proj=longlat +datum=WGS84") # WGS84
crs_proj <- CRS("+proj=utm +zone=17 +datum=WGS84") # UTM17

study_area_proj <- st_transform(study_area, crs_proj) # transform WGS84 -> UTM17 (unit m)
bbox <- st_bbox(study_area_proj)

r <- raster(xmn = bbox[1], xmx = bbox[3], ymn = bbox[2], ymx = bbox[4], crs = crs_proj, resolution = 500)
values(r) <- 1:ncell(r)

r <- terra::mask(r, study_area_proj) # Mask r by study area

plot(r) # plot r

r_idx <- values(r)[!is.na(values(r))] # cell indices where the study area occupies

# Generate raster centroids
r_centorid <- rasterToPoints(r, spatial = TRUE) %>% 
  st_as_sf() # 'layer' column = r_idx
colnames(r_centorid)[1] <- "cell_id"

# Part 1. Create predictor variables
## 1. Airport
airport <- read_sf('airport.shp') %>% st_transform(crs_proj)

dist_airport <- st_distance(airport, r_centorid) %>% as.data.frame() # distance between airport and centroid
dist_airport.df <- dist_airport %>% apply(MARGIN = 2, FUN = min) %>% as.data.frame() %>%
  tidyr::pivot_longer(cols = everything()) # select the nearest airport distance

dist_airport.df$name <- r_idx # change the name column as index

## raster map
dist_airport.r <- r # airport distance raster
values(dist_airport.r)[r_idx] <- dist_airport.df$value
plot(dist_airport.r)

pred_var.df <- dist_airport.df # pred_var.df is a final dataframe containing all predcitor variables
colnames(pred_var.df) <- c("cell_id", "dist_airport")

## 2. Emission of Title V
titleV <- read_sf('TitleV_emission_source.shp') %>% st_transform(crs_proj)

## [tons / km] (PM2.5 emission / distance)
dist_titleV <- st_distance(r_centorid, titleV) %>% as.data.frame() # each column = a facility 

dist_titleV_km <- dist_titleV %>% apply(MARGIN = 2, FUN = as.numeric) # convert numeric 
dist_titleV_km <- dist_titleV_km / 1000 # convert unit (m -> km)

pm_titleV <- titleV$pm2_5 %>% t() # empty vector for PM emission
emit_titleV <- dist_titleV_km # empty vecto for PM2.5 emmission / km

for (i in c(1:ncol(pm_titleV))) { # calculate [tons/km]
  for (j in c(1:nrow(dist_titleV_km))) {
    emit_titleV[j,i] <- pm_titleV[1,i] / dist_titleV_km[j,i]
  }
}

## sum up [tons/km] from all facilities
emit_titleV.df <- emit_titleV %>% apply(MARGIN = 1, FUN = sum) %>% 
  as.data.frame() %>% tidyr::pivot_longer(cols = everything())

emit_titleV.df$name <- r_idx # name = r_idx
colnames(emit_titleV.df) <- c("cell_id", "emission_titleV")

## raster map for Emission of TitleV (optional)
titleV.r <- r 
values(titleV.r)[r_idx] <- emit_titleV.df$emission_titleV
plot(titleV.r)
plot(titleV$geometry, add = TRUE)

pred_var.df <- left_join(pred_var.df, emit_titleV.df)

## 3. Number of DEC sources in each cell
dec_source <- read_sf('nys_sources_afc.shp') %>% st_transform(crs_proj)
dec_source_count <- cellFromXY(r, st_coordinates(dec_source)) %>% table() # how many DEC source are in cells in r?
dec_source.r <- r # raster for DEC variable
values(dec_source.r) <- rep(NA, ncell(dec_source.r)) # At first, fill NA values
values(dec_source.r)[as.numeric(names(dec_source_count))] <- dec_source_count # fill the cell value
plot(dec_source.r) # plot

hist(values(dec_source.r)) # histogram of the variable
dec_source.df <- data.frame(DEC_source = values(dec_source.r))
dec_source.df$DEC_source[is.na(dec_source.df$DEC_source)] <- 0

pred_var.df <- left_join(pred_var.df, emit_titleV.df) # join to pred_var.df

## 4. Road and highway
road <- read_sf("NYS_Streets_Intersect.shp") %>% st_transform(crs_proj)
highway <- read_sf("NYS_Streets_highway_Intersect.shp") %>% st_transform(crs_proj)

## Rasterization 10x10m
raster_template <- rast(ext(r), resolution = 10,
                       crs = st_crs(road)$wkt)

road.r <- rasterize(road, raster_template, touches = TRUE)
highway.r <- rasterize(highway, raster_template, touches = TRUE)

## Zonal statistic - Calculate percentage of road area in cells of r
zone.r <- r %>% as("SpatRaster") %>% terra::disagg(fact = 50) # change spatial resolution (disaggreating r (500x500) into 10x10)
road_percent <- terra::zonal(c(road.r, highway.r), zone.r, fun = "sum", na.rm = TRUE) # Sum cell values
colnames(road_percent) <- c("cell_id", "road_count", "highway_count")
road_percent$road_count[is.nan(road_percent$road_count)] <- 0 # Replace NaN value with 0
road_percent$highway_count[is.nan(road_percent$highway_count)] <- 0 # Replace NaN value with 0
road_percent$road_percent <- (road_percent$road_count*10*10) / (500*500) # percentage of occupied road area in 500x500m cell
road_percent$highway_percent <- (road_percent$highway_count*10*10) / (500*500) # percentage of occupied highway area in 500x500m cell

road_percent.r <- r # raster map for road
values(road_percent.r)[r_idx] <- road_percent$road_percent

highway_percent.r <- r # raster map for highway
values(highway_percent.r)[r_idx] <- road_percent$highway_percent

road_percent_selected <- road_percent %>% dplyr::select("cell_id", "road_percent", "highway_percent")
pred_var.df <- left_join(pred_var.df, road_percent_selected) # add to pred_var.df

## 5. Distance to railroad
railroad <- read_sf("Railroad.shp") %>% st_transform(crs_proj)
dist_railroad <- st_distance(r_centorid, railroad) %>% as.data.frame()

dist_railroad.df <- dist_railroad %>% apply(MARGIN = 1, FUN = min) %>% as.data.frame() %>%
  tidyr::pivot_longer(cols = everything())

dist_railroad.df$name <- r_idx
colnames(dist_railroad.df) <- c("cell_id", "dist_railroad")

railroad.r <- r # raster map for dist_railroad
values(railroad.r)[r_idx] <- dist_railroad.df$dist_railroad

pred_var.df <- left_join(pred_var.df, dist_railroad.df) # add to pred_var.df

## 6. Distance to Highway
dist_highway <- st_distance(r_centorid, highway) %>% as.data.frame()

dist_highway.df <- dist_highway %>% apply(MARGIN = 1, FUN = min) %>% as.data.frame() %>%
  tidyr::pivot_longer(cols = everything())

dist_highway.df$name <- r_idx
colnames(dist_highway.df) <- c("cell_id", "dist_highway")

pred_var.df <- left_join(pred_var.df, dist_highway.df) # add to pred_var.df


## 7. AADT
aadt <- read_sf("AADT_Intersect.shp") %>% st_transform(crs_proj) # AADT
truck_per <- read_sf("Percentage_Truck_Intersect.shp") %>% st_transform(crs_proj) #truck AADT

## Rasterize AADT and truckAADT by 10x10m 
aadt.r <- rasterize(aadt, raster_template, field = "AADT", fun = "mean", touches = TRUE)
truck_aadt.r <- rasterize(aadt, raster_template, field = "TruckAADT", fun = "mean", touches = TRUE)
truck_per.r <- rasterize(truck_per, raster_template,  field = "Percentile", fun = "mean", touches = TRUE) # truck percentage per segment

## Zonal statistics
aadt <- terra::zonal(c(aadt.r, truck_aadt.r, truck_per.r), zone.r, fun = "mean", na.rm = TRUE) # mean cell values
colnames(aadt) <- c("cell_id", "AADT", "truckAADT", "truck_percentile")
aadt$AADT[is.nan(aadt$AADT)] <- 0 # Replace NaN value with 0
aadt$truckAADT[is.nan(aadt$truckAADT)] <- 0 # Replace NaN value with 0
aadt$truck_percentile[is.nan(aadt$truck_percentile)] <- 0 # Replace NaN value with 0

pred_var.df <- left_join(pred_var.df, aadt) # add to pred_var.df

## 8. Land use
lc <- rast("C:/Users/jlee367/OneDrive - University at Buffalo/Project/LUR/Data/raster/lc_agg.tif")
r_cent_buf_250 <- st_buffer(r_centorid, dist = 250) # Buffer 250m

lc_buf <- terra::extract(lc, r_cent_buf_250) 
lc_buf.df <- lc_buf %>% group_by(ID, Legend) %>% summarise(count = n())
lc_buf.df$percent_area <- (lc_buf.df$count*900*100) / (250*250*pi) # percentage of area in 250m buffer

temp.df <- lc_buf.df %>% pivot_wider(names_from = Legend, values_from = percent_area) # table describing how many areas per landuse type
temp.df2 <- temp.df %>% dplyr::select(-count) %>% group_by(ID) %>% # group by (buffer_id = cell_id)
  summarise(across(everything(), \(x) sum(x, na.rm = TRUE)))

pred_var.df <- cbind(pred_var.df, temp.df2[2:13]) # add to pred_var.df

## 9. Canopy
canopy <- rast("C:/Users/jlee367/OneDrive - University at Buffalo/Project/LUR/Data/raster/nlcd_tcc_conus_2021_v2021-4_Ta3WVsmoT75aVYq6JVNf.tiff")
canopy[canopy > 100] <- NA # only 0~100 is valid

#### Buffer and average
canopy_buf <- terra::extract(canopy, r_cent_buf_250, fun = "mean", na.rm = TRUE) # extract canopy mean value in each 250m buffer
canopy <- canopy_buf$Layer1
pred_var.df <- cbind(pred_var.df, canopy) # add to pred_var.df


## 10. Distance to park and forest
green <- read_sf("park_forest.shp") %>% st_transform(crs_proj)
dist_green <- st_distance(r_centorid, green) %>% as.data.frame()
dist_green.df <- dist_green %>% apply(MARGIN = 1, FUN = min) %>% as.data.frame() %>%
  tidyr::pivot_longer(cols = everything())

dist_green.df$name <- r_idx
colnames(dist_green.df) <- c("cell_id", "dist_green")

pred_var.df <- left_join(pred_var.df, dist_green.df)


## 11. NDVI (0425-0503)
ndvi <- rast("C:/Users/jlee367/OneDrive - University at Buffalo/Project/LUR/Data/raster/NDVI_studyarea.tif")

ndvi_buf <- terra::extract(ndvi, r_cent_buf_250)
ndvi_buf.df <- ndvi_buf %>% group_by(ID) %>% summarise(NDVI = mean(NDVI_studyarea))
pred_var.df <- cbind(pred_var.df, ndvi_buf.df[2])


## 12. Population
pop <- read_sf("block.shp") %>% st_transform(crs_proj)
pop$pop_den <- pop$tot_popE / (st_area(pop)/1000) # people/km2

raster_template <- rast(ext(r), resolution = 100,
                        crs = st_crs(road)$wkt) # generate population raster by 100x100m

pop.r <- rasterize(pop, raster_template, field = "pop_den", touches = TRUE)

### Zonal statistic
zone.r <- r %>% as("SpatRaster") %>% terra::disagg(fact = 5) # change spatial resolution
pop_den <- terra::zonal(pop.r, zone.r, fun = "mean", na.rm = TRUE) # mean cell values
colnames(pop_den) <- c("cell_id", "pop_den")

pred_var.df <- cbind(pred_var.df, pop_den[2])

pop_den.r <- r # raster map
values(pop_den.r)[r_idx] <- pop_den$pop_den
plot(pop_den.r)

## 13. Building density
building <- read_sf("Building_Footprints.shp") %>% st_transform(crs_proj)
r_polygon <- rasterToPolygons(r) %>% st_as_sf() # This time, I converted the raster to polygon because building polygons are too coarse to convert it to grids
building_r_intersect <- st_intersection(r_polygon, building)
building_r_intersect$area <- st_area(building_r_intersect)
building_r_intersect.df <- building_r_intersect %>% as.data.frame() %>% group_by(layer) %>% 
  summarise(building_percent = as.numeric((sum(area)*100) / (500*500)))
colnames(building_r_intersect.df)[1] <- "cell_id"

pred_var.df <- left_join(pred_var.df, building_r_intersect.df)
pred_var.df$building_percent[is.na(pred_var.df$building_percent)] <- 0

## 14. DEM
dem <- rast("DEM_studyarea.tif")
dem_buf <- terra::extract(dem, r_cent_buf_250) 
dem_buf.df <- dem_buf %>% group_by(ID) %>% summarise(DEM = mean(DEM_studyarea))
pred_var.df <- cbind(pred_var.df, dem_buf.df[2])

## 15. Meteorological data
### Temperature
temp <- rast("air.sfc.2024.nc")
dim(temp) # 277 349 335

temp_agg <- mean(temp, na.rm = TRUE) - 273.15 # aggregate values by averaging in z-axis and convert unit Kelvin to Celcius

study_area2 <- st_transform(study_area, crs(temp)) # transform study area to temperature data

temp_agg_mask <- mask(temp_agg, study_area2)  # mask temperature data
temp_agg_extract <- st_transform(r_cent_buf_250, crs(temp_agg_mask)) %>% 
  terra::extract(temp_agg_mask, ., fun = "mean", na.rm = TRUE) # extract it by 250m buffer around cell centroids
colnames(temp_agg_extract)[2] <- "temp"

temp.r <- r # raster map for temperature 
values(temp.r)[r_idx] <- temp_agg_extract$temp
plot(temp.r) # plot

pred_var.df <- cbind(pred_var.df, temp_agg_extract[2])

### Pressure
pressure <- rast('pres.sfc.2024.nc')
dim(pressure) # 277 349 335

pressure_agg <- mean(pressure, na.rm = TRUE) # aggregate values in z-axis
pressure_agg_mask <- mask(pressure_agg, study_area2) # mask temperature data
pressure_agg_extract <- st_transform(r_cent_buf_250, crs(pressure_agg_mask)) %>% 
  terra::extract(pressure_agg_mask, ., fun = "mean", na.rm = TRUE) # extract it by 250m buffer around cell centroids

colnames(pressure_agg_extract)[2] <- "pressure"
pred_var.df <- cbind(pred_var.df, pressure_agg_extract[2])

### Relative humidity
humidity <- rast('rhum.2m.2024.nc')
dim(humidity) # 277 349 335
humidity_agg <- mean(humidity, na.rm = TRUE) # aggregate values in z-axis
humidity_agg_mask <- mask(humidity_agg, study_area2) # mask temperature data
humidity_agg_extract <- st_transform(r_cent_buf_250, crs(humidity_agg_mask)) %>% 
  terra::extract(humidity_agg_mask, ., fun = "mean", na.rm = TRUE) # extract it by 250m buffer around cell centroids

colnames(humidity_agg_extract)[2] <- "humidity"
pred_var.df <- cbind(pred_var.df, humidity_agg_extract[2])

### Boundary Mixing Layer
mixing_lyr <- rast('bmixl.hl1.2024.nc')
dim(mixing_lyr) # 277 349 335
mixing_lyr_agg <- mean(mixing_lyr, na.rm = TRUE) # aggregate values in z-axis
mixing_lyr_agg_mask <- mask(mixing_lyr_agg, study_area2) # crop and mask temperature data
mixing_lyr_agg_extract <- st_transform(r_cent_buf_250, crs(mixing_lyr_agg_mask)) %>% 
  terra::extract(mixing_lyr_agg_mask, ., fun = "mean", na.rm = TRUE) # extract it by 250m buffer around cell centroids

colnames(mixing_lyr_agg_extract)[2] <- "mixingLayer"

pred_var.df <- cbind(pred_var.df, mixing_lyr_agg_extract[2])

## Wind speed (monthly)
windspeed <- rast('wspd.10m.mon.ltm.nc')
dim(windspeed) # 277 349  12
windspeed_agg <- mean(windspeed, na.rm = TRUE) # aggregate values in z-axis
windspeed_agg_mask <- mask(windspeed_agg, study_area2) # crop and mask temperature data
windspeed_agg_extract <- st_transform(r_cent_buf_250, crs(windspeed_agg_mask)) %>% 
  terra::extract(windspeed_agg_mask, ., fun = "mean", na.rm = TRUE) # extract it by 250m buffer around cell centroids

colnames(windspeed_agg_extract)[2] <- "windspeed"
pred_var.df <- cbind(pred_var.df, windspeed_agg_extract[2])

## V-wind
vwind <- rast('vwnd.10m.2024.nc')
dim(vwind) # 277 349 335
vwind_agg <- mean(vwind, na.rm = TRUE) # aggregate values in z-axis
vwind_agg_mask <- mask(vwind_agg, study_area2) # crop and mask temperature data
vwind_agg_extract <- st_transform(r_cent_buf_250, crs(vwind_agg_mask)) %>% 
  terra::extract(vwind_agg_mask, ., fun = "mean", na.rm = TRUE) # extract it by 250m buffer around cell centroids

colnames(vwind_agg_extract)[2] <- "vwind"
pred_var.df <- cbind(pred_var.df, vwind_agg_extract[2])

## U-wind
uwind <- rast('uwnd.10m.2024.nc')
dim(uwind) # 277 349 335
uwind_agg <- mean(uwind, na.rm = TRUE) # aggregate values in z-axis
uwind_agg_mask <- mask(uwind_agg, study_area2) # crop and mask temperature data
uwind_agg_extract <- st_transform(r_cent_buf_250, crs(uwind_agg_mask)) %>% 
  terra::extract(uwind_agg_mask, ., fun = "mean", na.rm = TRUE) # extract it by 250m buffer around cell centroids

colnames(uwind_agg_extract)[2] <- "uwind"
pred_var.df <- cbind(pred_var.df, uwind_agg_extract[2])


# Part 2. Append Mobile monitoring data (2022-07-01 00:00:00 UTC ~ 2023-06-30 23:59:59 UTC)
mm_data <- read.csv("aclima_pro_buffalo_ambient_latest.csv") %>% 
  filter(modality ==  "pm_2.5") %>%  filter(metric ==  "median") %>% 
  dplyr::select(segment_id, value, geometry) %>% st_as_sf(wkt = "geometry") %>% 
  st_set_crs(4326) %>% st_transform(crs_proj)

mm_data.r <- rasterize(mm_data, r, field = "value", fun = mean, touches = TRUE)

projectRaster(mm_data.r, crs = crs_wgs84) %>% mapview() # plot 

mm_data.df <- mm_data.r %>% as.data.frame()
colnames(mm_data.df) <- "mobile_PM2.5"
mm_data.df$cell_id <- rownames(mm_data.df) %>% as.numeric()
mm_data.df <- mm_data.df[,c(2,1)]

final.df <- left_join(mm_data.df, pred_var.df)
final.df <- final.df[2:ncol(final.df)]

## Edit column names (e.g. deleting blank spaces and commas)
colnames(final.df) <- sub(" ", "_", colnames(final.df))
colnames(final.df) <- sub(",", "", colnames(final.df))
colnames(final.df) <- gsub("[[:space:]]+", "_", colnames(final.df))
colnames(final.df) <- gsub("/", "_", colnames(final.df))

final_mm.df <- final.df[!is.na(final.df$mobile_PM2.5),] # only leave cells having mobile monitoring data

write.csv(final_mm.df, "mobile_predictors_0116.csv", row.names = FALSE)

# final output: final_mm.df #


