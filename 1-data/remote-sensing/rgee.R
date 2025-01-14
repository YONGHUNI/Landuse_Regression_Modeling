library(reticulate)
#reticulate::miniconda_uninstall()

library(rgee)

#ee_install()
#conda_install( envname = "rgee", "earthengine-api")

#ee_install_upgrade()

#https://github.com/r-spatial/rgee/pull/364
#remotes::install_github(repo="r-spatial/rgee", ref = remotes::github_pull(364))




ee_check()


ee_Initialize(drive=TRUE, quiet = F)


library(sf)
library(tidyverse)
library(terra)

#Sys.setlocale('LC_ALL', 'ko_KR.UTF-8')

mask <- read_sf("./data/features/study_area.gpkg") %>% 
    st_make_valid() %>%
    st_transform(st_crs(4326)) %>%
    sf_as_ee()

region <- mask$geometry()$bounds()

region$aside(ee_print)


### MODIS NDVI

#https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD13Q1#description
dataset <- ee$ImageCollection('MODIS/061/MOD13Q1')$filterDate('2024-07-01', '2025-01-02')

buffalo_dataset <- dataset$filterBounds(mask)
    
buffalo_nvdi <- buffalo_dataset$select('NDVI')

buffalo_nvdi$aside(ee_print)





Map$centerObject(mask, zoom = 8)

Map$addLayer(
    eeObject = buffalo_nvdi$first()$clip(mask)
)




buffalo_nvdi$reproject("EPSG:4326",NULL,250)$clip(mask)$aside(ee_print)


# for retrieving a single image
ext_aoi <- ee_as_sf(mask) %>% st_transform(st_crs(0000))

ee_as_rast(
    buffalo_nvdi$first()$reproject("EPSG:4326",NULL,250),
    region = region,
    #dsn = NULL,
    via = "drive",
    #container = "rgee_backup",
    #scale = 232,
    #maxPixels = 11612160001,
    #lazy = FALSE,
    #public = TRUE,
    #add_metadata = TRUE,
    #timePrefix = TRUE,
    #quiet = FALSE,
    #...
) -> test







crop_ndvi <- terra::project(test,crs("epsg:XXXX"),threads = T) %>%
    terra::crop(ext(ext_aoi), snap="out") %>%
    terra::mask(mask = ext_aoi)



library(tmap)
library(stars)
tmap_mode("plot")

qtm(test)


tm_shape(crop_ndvi) +
    tm_raster(col_alpha = .45) +
    tm_shape(ext_aoi) +
    tm_polygons()


#scale factor 0.0001
terra::writeRaster(crop_ndvi*0.0001,"./data/MODIS+061+MOD13Q1+2018_01_01-2018_01_17.tiff",overwrite=T)




# for retrieving the entire image collection
buffalo_nvdi_wgs84 <- buffalo_nvdi$map(function(img){return(img$reproject("EPSG:4326",NULL,250))})


ee_imagecollection_to_local(
    ic = buffalo_nvdi_wgs84,
    region = region,
    dsn = "./data/rasters/NDVI/",
    via = "drive",
    container = "rgee_backup",
    scale = NULL,
    maxPixels = 1e+09,
    lazy = FALSE,
    public = TRUE,
    add_metadata = TRUE,
    timePrefix = TRUE,
    quiet = FALSE
)




### SRTM DEM

buffalo_dem <- ee$Image('CGIAR/SRTM90_V4')$clip(mask)



buffalo_dem$aside(ee_print)





Map$centerObject(mask, zoom = 8)

Map$addLayer(
    eeObject = buffalo_dem
)


ee_as_rast(
    buffalo_dem,
    region = region,
    #dsn = NULL,
    via = "drive",
    #container = "rgee_backup",
    #scale = 232,
    #maxPixels = 11612160001,
    #lazy = FALSE,
    #public = TRUE,
    #add_metadata = TRUE,
    #timePrefix = TRUE,
    #quiet = FALSE,
    #...
) -> rast_dem


writeRaster(rast_dem,"./data/rasters/DEM/SRTM/SRTM90_V4.tiff")




