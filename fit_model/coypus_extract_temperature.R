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

Dir.Data <- file.path(getwd(), "data/temp") # folder path for data

# get temp in January of every year
years <- c("2016", "2017", "2018", "2019", "2020", 
           "2021", "2022", "2023", "2024")

for (i in 1:length(years)) {
  toccitanie <- CDownloadS(
    Variable = "2m_temperature",
    DataSet = "reanalysis-era5-land-monthly-means",
    Type = "monthly_averaged_reanalysis",
    DateStart = paste0(years[i], "-01-01 00:00"),
    DateStop = paste0(years[i], "-01-31 23:00"),
    TZone = "Europe/Paris",
    TResolution = "month",
    TStep = 1,
    Extent = occitanie_raster, # our data.frame with Lat and Lon columns
    Dir = Dir.Data,
    FileName = paste0("Toccitanie_", years[i]),
    API_User = "olivier.gimenez@cefe.cnrs.fr",
    API_Key = "ccb31c25-7603-4cd7-8e88-97c8eb6e9cbd"
  )
  
  
  terra::writeRaster(x = toccitanie,
                     filename = paste0("data/temp/temp_", years[i], ".tif"),
                     overwrite = TRUE)
}

# toccitanie <- terra::rast("data/shp/temp.tif")

Plot.SpatRast(toccitanie_2016$`2m_temperature`) + 
  geom_sf(data = dpts_occitanie %>% 
            st_union() %>%
            st_transform(crs = 4326), 
             fill = NA)

# join all years together


# extract temperature per commune

temp_communes <- terra::extract(terra::rast("data/temp/temp_2016.tif"), 
                                dpts_occitanie_sub) %>%
  group_by(ID) %>%
  summarise(t2016 = weathermetrics::kelvin.to.celsius(
    mean(`2m_temperature`, na.rm = T)
    )
  )

for (i in 2:length(years)) {
  
  temp_df <- terra::extract(terra::rast(paste0("data/temp/temp_", 
                                               years[i], ".tif")), 
                            dpts_occitanie_sub) %>%
    group_by(ID) %>%
    summarise(!!paste0("t", years[i]) := weathermetrics::kelvin.to.celsius(
      mean(`2m_temperature`, na.rm = T)
      )
    )
  
  temp_communes <- left_join(temp_communes, temp_df)
  
}

temp_communes <- temp_communes %>% 
  mutate(commune = communes_clean)

saveRDS(temp_communes, "data/temp/temperature.rds")
