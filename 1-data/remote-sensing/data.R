library(tidyverse)
library(sf)
library(terra)

read_sf("./data/features/study_area.gpkg")

rast("./data/rasters/DEM/SRTM/SRTM90_V4.tiff") %>%
    terrain(v="slope", neighbors=8, unit="degrees") %>%
    writeRaster(.,"./data/rasters/DEM/SRTM/derivatives/aspect.tiff")

rast("./data/rasters/DEM/SRTM/SRTM90_V4.tiff") %>%
    terrain(v="aspect", neighbors=8, unit="degrees") %>%
    writeRaster(.,"./data/rasters/DEM/SRTM/derivatives/aspect.tiff")