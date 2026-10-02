#### Update Carbon Emissions Model and test results
### Started 24 October 2025 by Cat

### Looking at changes over time using data from 2008, 2014, 2020, and 2025

## 10 Sept 2026 by Cat
### UPDATED WITH LAI TO 2008, 2012, 2016, 2020, 2023, and 2015 for Holly Shelter instead of 2016

### housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## Load Libraries
library(dplyr)
library(ggplot2)
library(tidyterra)
library(tidyr)
library(sf)
library(terra)
library(FedData)
library(daymetr)
library(patchwork)
library(broom.mixed)
library(gridExtra)
library(exactextractr)
library(scPDSI)
library(SPEI)


### Set working directory
setwd("~/Documents/git/tnc_cprg/analyses/peatlands/")

################################################################################
## Load in full dataset to get z values
cleandat <- read.csv("output/clean_swails.csv")

#### Load in Shapefile
sf <- st_read("../../data/North Carolina/NC_Pocosin_Restoration_Sites_2026/") %>%
  mutate(Proj_Name = ifelse(Label == "Hydrologic Restoration Area 3", Label, Proj_Name))

### Source the Function
source("source/updated_peatprefeasibility.R")

dem_lidar <- rast("lidar_data/NC_DEM_10m.tif")

#create and save functions - key with FedData
not_all_na <- function(x) any(!is.na(x)) #for data cleaning

save <- sf %>% st_transform(4326)

# Define the sites with their names, filenames, and filter criteria
sites <- list(
  #list(name = "Van Swamp", filename = "vanswamp", proj_name = "Van Swamp Restoration"),
  #list(name = "Hofmann Forest", filename = "hofmannforest", proj_name = "Hofmann Forest Pocosin Rewetting"),
  #list(name = "Holly Shelter", filename = "hollyshelter", proj_name = "Holly Shelter Pocosin"),
  #list(name = "Angola Bay", filename = "angolabay", proj_name = "Angola Bay Restoration Area"),
  list(name = "PLNWR Area 3", filename = "plnwr_3", proj_name = "Hydrologic Restoration Area 3")
)


#foo <- sf %>%
#  filter(Proj_Name %in%  c("Van Swamp Restoration","Hofmann Forest Pocosin Rewetting",
#                           "Holly Shelter Pocosin", "Angola Bay Restoration Area"))

#st_centroid(foo$geometry)


# Define the years
years <- c(2008, 2012, 2015, 2016, 2020, 2023)  

# Loop through each site
for (site in sites) {
  # Filter the shapefile for this site
  sf_site <- sf %>%
    filter(Proj_Name == site$proj_name)
  
  title <- site$name
  filename_base <- site$filename
  standname <- names(sf_site)[1]
  
  # Loop through each year
  for (inputyear in years) {
    filename <- paste0(filename_base, "_", inputyear)
    
    ### SSURGO processing
    sfoutput <- peatpref(sf_site, title, filename, standname)
    
    ### Clean up the dataframe
    siteinformation <- sfoutput[[2]]
    
    drainedsite <- siteinformation %>% 
      ungroup() %>%
      mutate(tx = 0) %>%
      rename(wtd = watertable,
             vegtype = vegetation) %>%
      mutate(#wtd = max(wtd, na.rm = TRUE))
             wtd = weighted.mean(wtd[wtd>0], numacres[wtd>0]),
             wtd = wtd + 10)
    
    #### Get Climate data for potential restoration site
    sf_buffered <- sf_site %>%
      st_transform(32617) %>%
      st_buffer(50000) %>%
      st_transform(4326)
    
    climate <- get_ghcn_daily(
      template = sf_buffered, 
      label = "climate",
      elements = c("tmax", "tmin", "prcp"), 
      standardize = TRUE,
      years = c(2008, 2012, 2014, 2015, 2016, 2020, 2023, 2025),
      force.redo = TRUE
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
      summarize(totalGDD = max(GDD, na.rm = TRUE)) %>%
      ## Match station to geometry
      left_join(climate$spatial %>% rename(STATION = ID)) %>%
      st_as_sf(sf_column_name = "geometry") %>%
      select(-NAME)
    
    
    #### Prepare for PDSI
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
      complete(YEAR = 2008:2025, MONTH = 1:12)
    
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
      select(-MONTH)
    
    pdsi_results <- st_as_sf(pdsi_results) %>%
      filter(!is.na(scPDSI)) %>%
      select(scPDSI, geometry, YEAR) %>%
      distinct()
    
    
    ### Merge GDD information with sf
    drainedsite <- drainedsite %>%
      st_transform(crs(gdd)) %>%
      st_join(gdd %>% filter(YEAR == inputyear) %>% select(-YEAR, -STATION), join = st_nearest_feature)
    
    ### Merge scPDSI information with sf
    drainedsite <- drainedsite %>%
      st_transform(crs(pdsi_results)) %>%
      st_join(pdsi_results %>% filter(YEAR == inputyear) %>% select(-YEAR), join = st_nearest_feature)
    
    ### Get Elevation data
    ### Use DEM to get elevation information for each site
    #dem <- get_ned(
    #  template = drainedsite,
    #  label = "DEM",
    #  force.redo = TRUE
    #)
    
    drainedsite <- drainedsite %>%
      st_transform(crs(dem_lidar)) %>%
      st_buffer(5) %>%
      mutate(elevation = exactextractr::exact_extract(dem_lidar, geom,
                                                      "weighted_mean", weights="area"))
    
    
    
    #### Clean up the spatial data to get distance to coast
    ## Load in Coastal shapefile
    coast <- st_read("../input/ne_10m_coastline.shp") %>%
      st_transform(crs(drainedsite)) %>%
      st_make_valid()
    
    # Only compute distance to the nearest coastline feature
    coast_cropped <- st_crop(coast, st_buffer(st_make_valid(st_union(drainedsite)), 500000))
    dist <- st_distance(coast_cropped, drainedsite)
    drainedsite$distkm <- as.vector(apply(dist, 2, min)) / 1000
    
    ### Add in NDVI and LAI
    ## Load in raster file from GEE
    ndvirast <- rast(paste0("~/OneDrive - The Nature Conservancy/Extra CPRG Data/NC NDVI/ndvi", inputyear, ".tif"))
    lairast <- rast(paste0("~/OneDrive - The Nature Conservancy/Extra CPRG Data/MODIS_LAI/modis_lai_", inputyear, ".tif"))
    
    ## Append to dataframe
    drainedsite <-  drainedsite %>%
      mutate(ndvi = exact_extract(ndvirast, geom, "weighted_mean", weights = "area"),
             lai = exact_extract(lairast, geom, "weighted_mean", weights = "area"))
    
    
    ### Create new dataset to predict C flux after restoration
    restoredsites <- drainedsite %>%
      select(-totalGDD, -scPDSI) %>%
      st_transform(crs(gdd)) %>%
      st_join(gdd %>% filter(YEAR == inputyear) %>% select(-YEAR, -STATION), join = st_nearest_feature) %>%
      st_transform(crs(pdsi_results)) %>%
      st_join(pdsi_results %>% filter(YEAR == inputyear) %>% select(-YEAR), join = st_nearest_feature) %>%
      mutate(wtd = 20, tx = 1) %>%
      st_drop_geometry() %>%
      rbind(drainedsite %>% st_drop_geometry()) %>%
      mutate(vegetation = ifelse(vegtype == "Tree Cover", 1, 0),
             location = title,
             site = title,
             uniqueid = paste0(title, 1),
             soilresp = NA, 
             vegtype = NULL)
    
    
    write.csv(restoredsites, paste0("output/clean_", filename, "_restored.csv"), row.names = FALSE)
    
    print(paste0("Completed ", title, " for year ", inputyear))
    
  } # end year loop
  
} # end site loop