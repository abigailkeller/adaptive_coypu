# get temperature through kriging

library(tidyverse)
library(tidyterra)
library(lubridate)
library(stringi)
library(sf)
library(KrigR) # https://www.erikkusch.com/courses/krigr/

# communes
communes <- c(
  "baillargues", "candillargues", "entre-vignes", "la grande-motte",
  "lansargues", "lunel", "lunel-viel", "marsillargues", "mauguio", "perols",
  "saint-geniès-des-mourgues", "saint-just", "saint-nazaire-de-pézan",
  "saint-vincent-de-barbeyrargues", "valergues"
)
communes_clean <- stri_trans_general(communes, "Latin-ASCII")

# départements de la région
dpts_occitanie <- st_read(
  "data/communes-de-la-region-occitanie/georef-france-commune-millesime.shp"
  ) 
dpts_occitanie$dpts_clean <- tolower(
  stri_trans_general(dpts_occitanie$Nom_Officie.9, "Latin-ASCII")
  )

# subset to smaller region
dpts_occitanie_sub <- dpts_occitanie[dpts_occitanie$dpts_clean %in% 
                                       communes_clean, ]

# write to shapefile
# st_write(dpts_occitanie_sub, "data/shp/communes.shp")

# contours de la région
loc_site <- dpts_occitanie_sub %>% 
  st_union() %>%
  st_transform(crs = 4326)

st_bbox(loc_site)

occitanie_raster <- st_bbox(c(xmin = -0.35, 
                              xmax = 5, 
                              ymax = 45.1, 
                              ymin = 42), 
                            crs = st_crs(4326)) %>%
  st_as_sfc() %>%
  terra::vect() %>%
  terra::rast()

Dir.Base <- getwd() # identifying the current directory
Dir.Data <- file.path(Dir.Base, "data") # folder path for data

toccitanie <- CDownloadS(
  Variable = "2m_temperature",
  DataSet = "reanalysis-era5-land-monthly-means",
  Type = "monthly_averaged_reanalysis",
  DateStart = "2015-01-12 00:00",
  DateStop = "2016-03-31 23:00",
  TZone = "Europe/Paris",
  TResolution = "month",
  TStep = 1,
  Extent = occitanie_raster, # our data.frame with Lat and Lon columns
  Dir = Dir.Data,
  FileName = "Toccitanie",
  API_User = "olivier.gimenez@cefe.cnrs.fr",
  API_Key = "ccb31c25-7603-4cd7-8e88-97c8eb6e9cbd"
)


terra::writeRaster(x = toccitanie,
                  filename = "data/shp/temp.tif",
                  overwrite = TRUE)

toccitanie <- terra::rast("data/shp/temp.tif")

Plot.SpatRast(toccitanie$`2m_temperature_1`) + 
  geom_sf(data = dpts_occitanie %>% 
            st_union() %>%
            st_transform(crs = 4326), 
             fill = NA)

# extract temperature per commune

temp_communes <- terra::extract(toccitanie, dpts_occitanie_sub) %>%
  group_by(ID) %>%
  summarise(tdec = mean(`2m_temperature_11`, na.rm = T),
            tjan = mean(`2m_temperature_1`, na.rm = T),
            tfeb = mean(`2m_temperature_2`, na.rm = T),
            tmar = mean(`2m_temperature_3`, na.rm = T),
  )



temp_communes <- temp_communes %>% 
  mutate(tdec = weathermetrics::kelvin.to.celsius(tdec),
         tjan = weathermetrics::kelvin.to.celsius(tjan),
         tfeb = weathermetrics::kelvin.to.celsius(tfeb),
         tmar = weathermetrics::kelvin.to.celsius(tmar),
         commune = communes_clean)

saveRDS(temp_communes, "data/temperature.rds")
