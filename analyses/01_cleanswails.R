### Attempt to rerun Swails et al 2022 to incorporate more predictors
## This code will compile all the data and build a new model
### Started 23 October 2025 by Cat

### Consider adding in time since restoration

### Notes from 2 December 2025
## Look at GDDs rather than max temp
## explore precip variability or total precip during growing season
## Look at height above sea-level instead of distance to coast
## Add in Armstrong and Richardson data

### housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

### Load Libraries
library(sf)
library(dplyr)
library(tidyr)
library(mice)
library(brms)
library(pollen)
library(FedData)
library(terra)
### Get PDSI data for each site
library(SPEI)
library(scPDSI)
library(ggplot2)

### Set working directory
setwd("~/Documents/git/tnc_cprg/analyses/peatlands/")

set.seed(1221)

alldata <- read.csv("input/swailsdata.csv") %>%
  select(-CO2_Rs_check) %>%
  mutate(Plot = ifelse(Plot %in% c("D16_A", "D16_B"), "D16", Plot))

names(alldata) <- c("ref", "location", "status", "plot",
                    "year", "month", "date", "soilresp", "soiltemp", "wtl", "latitude", "longitude")

claytonblocks <- st_read("input/Clayton Blocks/armstrong_claytonblocks.shp") %>%
  mutate(plot = ifelse(Name == "Project Area 2", "C13",
                       ifelse(Name == "Project Area 3", "C14",
                              ifelse(Name == "Reference Site", "D16", NA))),
         ref = "Armstrong", 
         location = "PLNWR",
         Latitude = st_coordinates(st_centroid(geometry))[,2],
         Longitude = st_coordinates(st_centroid(geometry))[,1]) %>%
  select(ref, location, plot, Latitude, Longitude) %>%
  filter(!is.na(plot)) %>%
  st_drop_geometry()

alldata <- left_join(alldata, claytonblocks)



alldata <- alldata %>%
  mutate(date = NULL,
         locationyear = paste0(location, year),
         ### Fix error in Database
         precip = ifelse(location == "GDS", 1364,
                         ifelse(location == "HF", 1137,
                                ifelse(locationyear %in% c("PLNWR2011", "PLNWR2012", "PLNWR2013"), 1364,
                                       ifelse(locationyear %in% c("PLNWR2016", "PLNWR2017"), 1743,
                                              ifelse(locationyear %in% c("TLRP2007", "TLRP2008", "TLRP2009"), 1151, 1352))))),
         maxt = ifelse(location == "GDS", 23,
                       ifelse(location == "HF", 23.6,
                              ifelse(locationyear %in% c("PLNWR2011", "PLNWR2012", "PLNWR2013"), 21.3,
                                     ifelse(locationyear %in% c("PLNWR2016", "PLNWR2017"), 22.6,
                                            ifelse(locationyear %in% c("TLRP2007", "TLRP2008", "TLRP2009"), 22.9, 21.8))))),
         mint = ifelse(location == "GDS", 11.1,
                       ifelse(location == "HF", 11.9,
                              ifelse(locationyear %in% c("PLNWR2011", "PLNWR2012", "PLNWR2013"), 10,
                                     ifelse(locationyear %in% c("PLNWR2016", "PLNWR2017"), 11.6,
                                            ifelse(locationyear %in% c("TLRP2007", "TLRP2008", "TLRP2009"), 10.9, 10.2))))),
         precip = ifelse(location == "GDS", 1364,
                         ifelse(location == "HF", 1137,
                                ifelse(locationyear %in% c("PLNWR2011", "PLNWR2012", "PLNWR2013"), 1364,
                                       ifelse(locationyear %in% c("PLNWR2016", "PLNWR2017"), 1743,
                                              ifelse(locationyear %in% c("TLRP2007", "TLRP2008", "TLRP2009"), 1151, 1352))))),
         tx = ifelse(status == "Drained", 0, 1),
         latitude = ifelse(location %in% c("TLRP"), 36.62, 
                           ifelse(location %in% c("HF"), 34.8250, 
                                  ifelse(ref == "Armstrong", Latitude, latitude))), 
         longitude = ifelse(location %in% c("TLRP"), -76.47,
                            ifelse(location == "HF", -77.3220, 
                                   ifelse(ref == "Armstrong", Longitude, longitude))),
         vegsite = paste0(location, status),
         vegetation = ifelse(vegsite %in% c("HFDrained", "PLNWRDrained"), 0, 1)) %>% ## 0 is herbaceous, 1 is tree cover
         select(-Latitude, -Longitude)

################################################################################
################################################################################
### Update the climate data using GHCN data instead
plots <- st_as_sf(alldata, 
                  coords = c("longitude", "latitude"), crs = 4326) 

climate <- get_ghcn_daily(
  template = plots, 
  label = "climate",
  elements = c("tmax", "tmin", "prcp"), 
  standardize = TRUE,
  years = c(2007:2017)
)

### TMAX
tmax <- as.data.frame(do.call(bind_rows, climate$tabular))$TMAX %>%
  pivot_longer(cols = c(D1:D31), names_to = "DAY", values_to = "TMAX") %>%
  mutate(TMAX = TMAX/10, ### Convert to degrees C
         DAY = as.numeric(gsub("D", "", DAY))) %>%
  filter(!is.na(TMAX)) 

### TMIN
tmin <- as.data.frame(do.call(bind_rows, climate$tabular))$TMIN %>%
  pivot_longer(cols = c(D1:D31), names_to = "DAY", values_to = "TMIN") %>%
  mutate(TMIN = TMIN/10, ### Convert to degrees C
         DAY = as.numeric(gsub("D", "", DAY))) %>%
  filter(!is.na(TMIN)) 

### Get GDDs
gdd <- left_join(tmax, tmin) %>%
  mutate(GDD = pmax(0, (TMAX + TMIN) / 2 - 5)) %>%
  group_by(STATION, YEAR,) %>%
  mutate(GDD = cumsum(GDD)) %>%
  group_by(STATION, YEAR, MONTH) %>%
  summarize(GDD = max(GDD)) %>%
  ungroup() %>%
  group_by(STATION, YEAR) %>%
  mutate(totalGDD = max(GDD, na.rm = TRUE))

## Get Growing Season GDDs
gsgdd <- left_join(tmax, tmin) %>%
  mutate(GDD = pmax(0, (TMAX + TMIN) / 2 - 5)) %>%
  group_by(STATION, YEAR,) %>%
  mutate(GDD = cumsum(GDD))

## Get Monthly Temps
meanmax <- tmax %>%
  group_by(STATION, YEAR, MONTH) %>%
  summarize(TMAX = mean(TMAX, na.rm = TRUE))

meanmin <- tmin %>%
  group_by(STATION, YEAR, MONTH) %>%
  summarize(TMIN = mean(TMIN, na.rm = TRUE))

### PRECIP
prcp <- as.data.frame(do.call(bind_rows, climate$tabular))$PRCP %>%
  pivot_longer(cols = c(D1:D31), names_to = "DAY", values_to = "PRCP") %>%
  mutate(PRCP = PRCP/10, ### Convert to mm 
         DAY = as.numeric(gsub("D", "", DAY))) %>%
  filter(!is.na(PRCP)) 

## Monthly Precip
totalprcp <- prcp %>%
  group_by(STATION, YEAR, MONTH) %>%
  summarize(PRCP = sum(PRCP, na.rm = TRUE))

### Growing Season Total Precip
gsprcp <- prcp %>%
  filter(MONTH %in% c("04", "05", "06",
                      "07", "08", "09", "10")) %>%
  group_by(STATION, YEAR) %>%
  summarize(PRCP = sum(PRCP, na.rm = TRUE))


### Put all climate data together
clim <- left_join(gdd, meanmax) %>%
  left_join(meanmin) %>%
  left_join(gsprcp) %>%
  ## Match station to geometry
  left_join(climate$spatial %>% rename(STATION = ID)) %>%
  ## Match formatting in swails dataset
  mutate(month = month.name[as.numeric(MONTH)],
         year = as.numeric(YEAR)) %>%
  ungroup() %>%
  mutate(geometry = st_transform(geometry, crs = 4326)) %>%
  select(year, month, totalGDD, TMAX, TMIN, PRCP, geometry) %>%
  mutate(month.year = paste(year, month))

clim <- st_as_sf(clim, sf_column_name = "geometry")


#### NOTE 15 April 2026 by Cat
## Need to create a buffer around the climate data and then intersect the Swails data points with
# the closest weather station


alldata2 <- alldata %>%
  mutate(month.year = paste(year, month)) %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(clim$geometry)) %>%
  st_transform(4326)


### Join the site data to the nearest weather station for the month and year
swails_climate <- do.call(rbind, lapply(unique(alldata2$month.year), function(cat) {
  st_join(
    subset(alldata2, month.year == cat), 
    subset(clim, month.year == cat), 
    join = st_nearest_feature
  )})) %>%
  select(-month.year.x, -month.y, -year.y, -month.year.y) %>%
  rename(month = month.x, year = year.x)

#st_write(swails_climate, "output/swails_climate.gpkg")


################################################################################
################################################################################
### Calculate PDSI for each site and month

# Calculate PET using latitude and SPEI package
pdsi <- totalprcp %>%
  left_join(tmax %>% group_by(STATION, YEAR, MONTH) %>% 
              summarize(TMAX = mean(TMAX, na.rm = TRUE))) %>%
  left_join(tmin %>% group_by(STATION, YEAR, MONTH) %>% 
              summarize(TMIN = mean(TMIN, na.rm = TRUE))) %>%
  ## Match station to geometry
  left_join(climate$spatial %>% rename(STATION = ID)) %>%
  mutate(latitude = st_coordinates(geometry)[,"Y"]) %>%
  group_by(STATION, YEAR, MONTH) %>%
  mutate(PET = hargreaves(TMIN, TMAX, lat = latitude, na.rm = TRUE),
         MONTH = as.numeric(MONTH),
         YEAR = as.numeric(YEAR)) %>%
  arrange(STATION, YEAR, MONTH)

### pdsi function needs a complete dataset to calculate for each station
# remove stations without some month's worth of data
pdsi <- pdsi %>%
  ungroup() %>%
  group_by(STATION) %>%
  complete(YEAR = 2007:2017, MONTH = 1:12)

run_scpdsi <- function(station_df) {
  
  # Get January start for each station
  start_year <- station_df$YEAR[1]
  
  # All months in dataset need to be complete
  P  <- station_df$PRCP   # precipitation (mm)
  PE <- station_df$PET    # pre-computed PET (mm)
  
  # Run self-calibrating PDSI
  result <- pdsi(P  = P,
                 PE = PE,
                 start = start_year,
                 sc = TRUE)   # sc = TRUE → scPDSI; FALSE → classic PDSI
  
  # pdsi() returns a list; $X is the PDSI time series
  n <- length(result$X)
  
  # Rebuild a year/month index that matches the output length
  all_months <- seq_len(n)
  years  <- start_year + (all_months - 1) %/% 12
  months <- ((all_months - 1) %% 12) + 1
  
  data.frame(
    STATION = station_df$STATION[1],
    YEAR    = years,
    MONTH   = months,
    scPDSI  = as.numeric(result$X),
    stringsAsFactors = FALSE
  )
}

pdsi_results <- pdsi %>%
  group_by(STATION) %>%
  # list of per-station data frames
  group_split() %>%  
  # run scPDSI for each
  purrr::map(run_scpdsi) %>%   
  # stack into one data frame
  bind_rows() %>%
  ## Match station to geometry
  left_join(climate$spatial %>% rename(STATION = ID)) %>%
  mutate(MONTH = month.name[as.numeric(MONTH)],
         month.year = paste(YEAR, MONTH)) %>%
  select(-MONTH, -YEAR)

pdsi_results <- st_as_sf(pdsi_results)

swails_climate$month.year <- paste(swails_climate$year, swails_climate$month)

swails_climate <- do.call(rbind, lapply(unique(swails_climate$month.year), function(cat) {
  st_join(
    subset(swails_climate, month.year == cat), 
    subset(pdsi_results, month.year == cat), 
    join = st_nearest_feature
  )})) %>%
  select(-month.year.x, -month.year.y, -STATION, -NAME) 



################################################################################
################################################################################
### Use Lidar data to get elevation information for each site

dem_lidar <- rast("lidar_data/NC_DEM_10m.tif")

swails_climate <- swails_climate %>%
  st_transform(crs(dem_lidar)) %>%
  st_buffer(5) %>%
  mutate(elevation = exactextractr::exact_extract(dem_lidar, geometry,
      "weighted_mean", weights = "area"))

#dem <- get_ned(
#  template = plots,
#  label = "DEM"
#)

#swails_climate <- swails_climate %>%
#  st_transform(crs(dem)) %>%
#  st_buffer(5) %>%
#  mutate(elevation = exactextractr::exact_extract(dem, geometry,
#                                                  "weighted_mean", weights="area"))


################################################################################
################################################################################
if(FALSE){ ## Not granular enough
  ### Use NLCD to get more refined vegetation information
  nlcd <- get_nlcd(
    template = plots,
    label = "NLCD",
    force.redo = TRUE
  )
  
  
  swails_climate <- swails_climate %>%
    st_transform(crs(nlcd)) %>%
    mutate(nlcd = exactextractr::exact_extract(nlcd, geometry,"mode"))
}

if(TRUE){ ## Try LAI instead
### Using NDVI instead
## Load in raster file from GEE
ndvi07 <- rast("~/OneDrive - The Nature Conservancy/Extra CPRG Data/NC NDVI/ndvi2007.tif")
ndvi08 <- rast("~/OneDrive - The Nature Conservancy/Extra CPRG Data/NC NDVI/ndvi2008.tif")
ndvi09 <- rast("~/OneDrive - The Nature Conservancy/Extra CPRG Data/NC NDVI/ndvi2009.tif")
ndvi11 <- rast("~/OneDrive - The Nature Conservancy/Extra CPRG Data/NC NDVI/ndvi2011.tif")
ndvi12 <- rast("~/OneDrive - The Nature Conservancy/Extra CPRG Data/NC NDVI/ndvi2012.tif")
ndvi13 <- rast("~/OneDrive - The Nature Conservancy/Extra CPRG Data/NC NDVI/ndvi2013.tif")
ndvi15 <- rast("~/OneDrive - The Nature Conservancy/Extra CPRG Data/NC NDVI/ndvi2015.tif")
ndvi16 <- rast("~/OneDrive - The Nature Conservancy/Extra CPRG Data/NC NDVI/ndvi2016.tif")
ndvi17 <- rast("~/OneDrive - The Nature Conservancy/Extra CPRG Data/NC NDVI/ndvi2017.tif")

## Append to dataframe
swails_climate <-  swails_climate %>%
  mutate(ndvi = case_when(
    year == 2007  ~ exactextractr::exact_extract(ndvi07, geometry, "weighted_mean", weights = "area"),
    year == 2008  ~ exactextractr::exact_extract(ndvi08, geometry, "weighted_mean", weights = "area"),
    year == 2009  ~ exactextractr::exact_extract(ndvi09, geometry, "weighted_mean", weights = "area"),
    year == 2011  ~ exactextractr::exact_extract(ndvi11, geometry, "weighted_mean", weights = "area"),
    year == 2012  ~ exactextractr::exact_extract(ndvi12, geometry, "weighted_mean", weights = "area"),
    year == 2013  ~ exactextractr::exact_extract(ndvi13, geometry, "weighted_mean", weights = "area"),
    year == 2015  ~ exactextractr::exact_extract(ndvi15, geometry, "weighted_mean", weights = "area"),
    year == 2016  ~ exactextractr::exact_extract(ndvi16, geometry, "weighted_mean", weights = "area"),
    year == 2017  ~ exactextractr::exact_extract(ndvi17, geometry, "weighted_mean", weights = "area"),
  ))
}


### Using LAI instead
## Load in raster file from GEE
lai07 <- rast("~/OneDrive - The Nature Conservancy/Extra CPRG Data/MODIS_LAI/modis_lai_2007.tif")
lai08 <- rast("~/OneDrive - The Nature Conservancy/Extra CPRG Data/MODIS_LAI/modis_lai_2008.tif")
lai09 <- rast("~/OneDrive - The Nature Conservancy/Extra CPRG Data/MODIS_LAI/modis_lai_2009.tif")
lai11 <- rast("~/OneDrive - The Nature Conservancy/Extra CPRG Data/MODIS_LAI/modis_lai_2011.tif")
lai12 <- rast("~/OneDrive - The Nature Conservancy/Extra CPRG Data/MODIS_LAI/modis_lai_2012.tif")
lai13 <- rast("~/OneDrive - The Nature Conservancy/Extra CPRG Data/MODIS_LAI/modis_lai_2013.tif")
lai15 <- rast("~/OneDrive - The Nature Conservancy/Extra CPRG Data/MODIS_LAI/modis_lai_2015.tif")
lai16 <- rast("~/OneDrive - The Nature Conservancy/Extra CPRG Data/MODIS_LAI/modis_lai_2016.tif")
lai17 <- rast("~/OneDrive - The Nature Conservancy/Extra CPRG Data/MODIS_LAI/modis_lai_2017.tif")

## Append to dataframe
swails_climate <-  swails_climate %>%
  mutate(lai = case_when(
    year == 2007  ~ exactextractr::exact_extract(lai07, geometry, "weighted_mean", weights = "area"),
    year == 2008  ~ exactextractr::exact_extract(lai08, geometry, "weighted_mean", weights = "area"),
    year == 2009  ~ exactextractr::exact_extract(lai09, geometry, "weighted_mean", weights = "area"),
    year == 2011  ~ exactextractr::exact_extract(lai11, geometry, "weighted_mean", weights = "area"),
    year == 2012  ~ exactextractr::exact_extract(lai12, geometry, "weighted_mean", weights = "area"),
    year == 2013  ~ exactextractr::exact_extract(lai13, geometry, "weighted_mean", weights = "area"),
    year == 2015  ~ exactextractr::exact_extract(lai15, geometry, "weighted_mean", weights = "area"),
    year == 2016  ~ exactextractr::exact_extract(lai16, geometry, "weighted_mean", weights = "area"),
    year == 2017  ~ exactextractr::exact_extract(lai17, geometry, "weighted_mean", weights = "area"),
  ))


################################################################################
################################################################################
if(TRUE){ ## Use both elevation and distance to coast
  #GDS: 36.61742381088771, -76.46975496072177
  #PLNWR: 35° 37′–35° 44′ N, 76° 27′–76° 35′ W; 35.91558356544023, -76.25380570307722
  #HF: 34.87657994186733, -77.36957038315614
  # TLRP: 36.61742381088771, -76.46975496072177
  
  coast <- st_read("../input/ne_10m_coastline.shp") %>%
    st_transform(crs(swails_climate))
  
  # Only compute distance to the nearest coastline feature
  coast_cropped <- st_crop(coast, st_buffer(st_union(swails_climate), 500000))
  dist <- st_distance(coast_cropped, swails_climate)
  swails_climate$distkm <- as.vector(apply(dist, 2, min)) / 1000
  
}

################################################################################
################################################################################
if(FALSE){### Not enough granular data
  
  ### Determine where there is soil temperature data for this region
  # Get the raw GHCN inventory for your area
  inventory <- get_ghcn_inventory(
    template = plots_buffer
  )
  
  # What do we have?
  table(sort(inventory$ELEMENT))
  
  plots_buffer <- plots %>%
    st_transform(32617) %>%   # project to meters (UTM zone 17N for NC)
    st_buffer(50000) %>%      # 50km buffer
    st_transform(4326)
  
  
  ### Update 15 April 2026 - GHCN soil temperatue values are not returning
  soiltemp <- get_ghcn_daily(
    template = plots_buffer, 
    label = "soils",
    ### Selected based on above outputs
    elements = c("SN52", "SX52"), 
    standardize = TRUE,
    years = c(2007:2017),
    force.redo = TRUE
  )
  
  ### Soil Temps
  soiltemp <- as.data.frame(do.call(bind_rows, soils$tabular))$SN11 %>%
    pivot_longer(cols = c(D1:D31), names_to = "DAY", values_to = "SMIN") %>%
    mutate(SMIN = SMIN/10, ### Convert to degrees C
           DAY = as.numeric(gsub("D", "", DAY))) %>%
    filter(!is.na(SMIN)) %>%
    left_join(
      as.data.frame(do.call(bind_rows, soils$tabular))$SX11 %>%
        pivot_longer(cols = c(D1:D31), names_to = "DAY", values_to = "SMIN") %>%
        mutate(SMIN = SMIN/10, ### Convert to degrees C
               DAY = as.numeric(gsub("D", "", DAY))) %>%
        filter(!is.na(SMIN))
    )
}

if(TRUE){### There is not enough GHCN data for the region
  ### Update 15 April 2026 - GHCN soil temperature values are not returning
  ### Impute missing soil temps
  ### Drop geometry before imputation — sf geometry column causes mice to fail
  swails_climate_df <- swails_climate %>% 
    st_drop_geometry() %>%
    mutate(wtl = ifelse(ref == "Helton", wtl * -1, wtl))
  imp <- mice(swails_climate_df, m = 1, print = FALSE, method = "norm.predict")
  
  soiltempmod <- brm_multiple(soiltemp ~ soilresp + totalGDD + scPDSI + elevation + tx + month, 
                              data = imp, chains = 2)
  
  save(soiltempmod, file = "models/soiltemp_impute.Rdata")
  
  swails_climate$soiltemp.clean <- predict(soiltempmod, newdat = swails_climate,)[,"Estimate"]
  
  plot(swails_climate$soiltemp.clean, swails_climate$soiltemp)
  ### Impute missing water table depths
  wtlmod <- brm_multiple(wtl ~ soiltemp + soilresp + totalGDD + scPDSI + elevation + lai + tx + month,
                         data = imp, chains = 2)
  
  swails_climate$wtl.clean <- predict(wtlmod, newdat = swails_climate)[,"Estimate"]
  
  
  
  ggplot(swails_climate %>% 
           mutate(needclean = ifelse(is.na(wtl), "imputed", "measured")
                  ), 
         aes(x = wtl.clean, y = wtl, fill = needclean, col = needclean)) + geom_point()
}

################################################################################
################################################################################

cleandat <- swails_climate %>%
  mutate(wtl = ifelse(is.na(wtl), wtl.clean, 
                      ifelse(ref == "Helton", wtl * -1, wtl)),
         #wtd.old = wtl * -1,
         #wtd = ifelse(ref == "Helton", wtl * -1, wtl),
         wtd = wtl * -1,
         wtl = NULL,
         wtl.clean = NULL,
         wtd.z = (wtd - mean(wtd, na.rm=TRUE)) / (2 * sd(wtd, na.rm=TRUE)),
         #wtd.old.z = (wtd.old - mean(wtd.old, na.rm=TRUE)) / (2 * sd(wtd.old, na.rm=TRUE)),
         soiltemp = ifelse(is.na(soiltemp), soiltemp.clean, soiltemp),
         soiltemp.clean = NULL,
         soiltemp.z = (soiltemp - mean(soiltemp, na.rm=TRUE)) / (2 * sd(soiltemp, na.rm=TRUE)),
         #precip.z = (precip - mean(precip, na.rm=TRUE)) / (2 * sd(precip, na.rm=TRUE)),
         pdsi.z = (scPDSI - mean(scPDSI, na.rm=TRUE)) / (2 * sd(scPDSI, na.rm=TRUE)),
         #maxt.z = (maxt - mean(maxt, na.rm=TRUE)) / (2 * sd(maxt, na.rm=TRUE)),
         #mint.z = (mint - mean(mint, na.rm=TRUE)) / (2 * sd(mint, na.rm=TRUE)),
         gdd.z = (totalGDD - mean(totalGDD, na.rm=TRUE)) / (2 * sd(totalGDD, na.rm=TRUE)),
         dist.z = (distkm - mean(distkm, na.rm=TRUE)) / (2 * sd(distkm, na.rm=TRUE)),
         elev.z = (elevation - mean(elevation, na.rm=TRUE)) / (2 * sd(elevation, na.rm=TRUE)),
         tx.z = (tx - mean(tx, na.rm=TRUE)) / (2 * sd(tx, na.rm=TRUE)),
         year.z = (year - mean(year, na.rm=TRUE)) / (2 * sd(year, na.rm=TRUE)),
         #veg.z = (vegetation - mean(vegetation, na.rm=TRUE)) / (2 * sd(vegetation, na.rm=TRUE))
         ndvi.z = (ndvi - mean(ndvi, na.rm = TRUE)) / (2 * sd(ndvi, na.rm = TRUE)),
         lai.z = (lai - mean(lai, na.rm = TRUE)) / (2 * sd(lai, na.rm = TRUE))
         ) %>%
  st_drop_geometry() %>%
  left_join(alldata %>% select(ref, location, status, plot, latitude, longitude) %>% distinct())

### Check how many NAs
foo <- cleandat[complete.cases(cleandat),] ## only 8, keeping full dataset for now

### Write the output
write.csv(cleandat, "output/clean_swails.csv", row.names = FALSE)




