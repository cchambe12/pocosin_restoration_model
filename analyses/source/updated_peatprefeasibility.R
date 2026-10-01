##### Build Peat Assessment function
### Started 24 October 2025 by Cat - UPDATED!

### Consider incorporating National Hyrdrography Dataset
#https://www.usgs.gov/national-hydrography/national-hydrography-dataset

### Source polygon splitter function
source("../North Carolina/source/splitpolygons.R")

################################################################################
##### Issue with FedData using outdated links for SSURGO. 
## Below are some workarounds

# Patch 1: inventory shapefile download
fixed_download_ssurgo_inventory <- function(raw.dir, ...) {
  url <- "https://websoilsurvey.nrcs.usda.gov/DataAvailability/SoilDataAvailabilityShapefile.zip"
  destdir <- raw.dir
  FedData:::download_data(url = url, destdir = destdir, ...)
  return(normalizePath(paste(destdir, "/SoilDataAvailabilityShapefile.zip", sep = "")))
}

# Patch 2: individual study area zip download (nc=FALSE to avoid caching corrupt zips)
fixed_download_ssurgo_study_area <- function(area, date, raw.dir) {
  url <- paste("https://websoilsurvey.nrcs.usda.gov/DSD/Download/Cache/SSA/wss_SSA_",
               area, "_[", date, "].zip", sep = "")
  destdir <- raw.dir
  FedData:::download_data(url = url, destdir = destdir, nc = FALSE)
  return(normalizePath(paste(destdir, "/wss_SSA_", area, "_[", date, "].zip", sep = "")))
}

# Patch 3: fix MULTISURFACE geometry parsing from WFS
fixed_get_ssurgo_inventory <- function(template = NULL, raw.dir) {
  if (!is.null(template)) {
    template %<>%
      FedData:::template_to_sf() %>%
      sf::st_transform(4326)
  }
  
  if (
    !is.null(template) &&
    httr::status_code(
      httr::RETRY(
        verb = "GET",
        url = "https://sdmdataaccess.nrcs.usda.gov/Spatial/SDMWGS84Geographic.wfs"
      )
    ) == 200
  ) {
    bounds <- template %>% sf::st_bbox() %>% sf::st_as_sfc()
    
    if ((sf::st_bbox(template)[["xmax"]] - sf::st_bbox(template)[["xmin"]]) > 1 |
        (sf::st_bbox(template)[["ymax"]] - sf::st_bbox(template)[["ymin"]]) > 1) {
      bounds %<>% sf::st_intersection(FedData:::grid)
    }
    
    SSURGOAreas <- bounds %>%
      purrr::map_dfr(function(x) {
        bound <- x %>% sf::st_bbox()
        if (identical(bound["xmin"], bound["xmax"])) bound["xmax"] <- bound["xmax"] + 1e-04
        if (identical(bound["ymin"], bound["ymax"])) bound["ymax"] <- bound["ymax"] + 1e-04
        bbox.text <- paste(bound, collapse = ",")
        temp.file <- paste0(tempdir(), "/soils.gml")
        
        httr::RETRY(
          verb = "GET",
          url = "https://sdmdataaccess.nrcs.usda.gov/Spatial/SDMWGS84Geographic.wfs",
          query = list(
            Service = "WFS", Version = "1.1.0", Request = "GetFeature",
            Typename = "SurveyAreaPoly", BBOX = bbox.text,
            SRSNAME = "EPSG:4326", OUTPUTFORMAT = "GML3"
          ),
          httr::write_disk(temp.file, overwrite = TRUE)
        )
        
        tryCatch(
          suppressMessages(suppressWarnings(
            sf::read_sf(temp.file, drivers = "GML", type = 3) %>%  # type=3 forces MULTIPOLYGON
              dplyr::mutate(saverest = as.Date(
                lubridate::parse_date_time(saverest, orders = "b d Y HMOp", locale = "en_US")
              )) %>%
              sf::st_drop_geometry()
          )),
          error = function(e) return(NULL)
        )
      }) %>%
      dplyr::distinct() %>%
      dplyr::arrange(areasymbol)
  } else {
    tmpdir <- tempfile()
    if (!dir.create(tmpdir)) stop("failed to create my temporary directory")
    file <- FedData:::download_ssurgo_inventory(raw.dir = raw.dir)
    utils::unzip(file, exdir = tmpdir)
    SSURGOAreas <- sf::read_sf(normalizePath(tmpdir), layer = "soilsa_a_nrcs")
    if (!is.null(template)) {
      SSURGOAreas %<>%
        sf::st_make_valid() %>%
        sf::st_intersection(sf::st_transform(template, sf::st_crs(SSURGOAreas)))
    }
    unlink(tmpdir, recursive = TRUE)
  }
  
  if (0 %in% SSURGOAreas$iscomplete) {
    warning("Some of the soil surveys in your area are unavailable.\n",
            paste0(as.vector(SSURGOAreas[SSURGOAreas$iscomplete == 0, ]$areasymbol), collapse = "\n"))
  }
  
  return(SSURGOAreas)
}


# Wrapper that patches, runs get_ssurgo, then restores originals (even on error)
get_ssurgo_fixed <- function(template, label, ...) {
  # Save originals
  orig_download_ssurgo_inventory  <- FedData:::download_ssurgo_inventory
  orig_download_ssurgo_study_area <- FedData:::download_ssurgo_study_area
  orig_get_ssurgo_inventory       <- FedData:::get_ssurgo_inventory
  
  # Apply patches
  assignInNamespace("download_ssurgo_inventory",  fixed_download_ssurgo_inventory,  ns = "FedData")
  assignInNamespace("download_ssurgo_study_area", fixed_download_ssurgo_study_area, ns = "FedData")
  assignInNamespace("get_ssurgo_inventory",       fixed_get_ssurgo_inventory,       ns = "FedData")
  
  # Restore originals when function exits (even on error)
  on.exit({
    assignInNamespace("download_ssurgo_inventory",  orig_download_ssurgo_inventory,  ns = "FedData")
    assignInNamespace("download_ssurgo_study_area", orig_download_ssurgo_study_area, ns = "FedData")
    assignInNamespace("get_ssurgo_inventory",       orig_get_ssurgo_inventory,       ns = "FedData")
  })
  
  get_ssurgo(template = template, label = label, ...)
}


### Begin function
peatpref <- function(sf, title, filename, standname){
  
  ### Save the original sf for mapping
  sfsave <- sf
  
  sfbuff <- sf %>%
    st_make_valid() %>%
    st_buffer(1000) %>%
    st_make_valid()
  
  ## Start by splitting up the shapefile
  totalacres <- sum(units::set_units(st_area(sfbuff), "acre")) #Take care of units
  
  min_map_unit = as.vector(totalacres/((0.5/100) * totalacres))
  
  if(nrow(sf) < min_map_unit){
    
    min_map_unit = as.numeric(min_map_unit)
    #Dissolve shapefile into single record
    totshp <- sfbuff %>% 
      group_by() %>% 
      summarise() %>%
      st_make_valid()
    
    shp_acres <- units::set_units(st_area(totshp), "acre") #Take care of units
    shp_segments <- round(shp_acres / min_map_unit) #Divide shapefile acreage by desired polygon subdivision acreage to get total number of polygons for subdivision
    
    totshp <- split_poly(totshp, shp_segments)
    totshp$AcresMade <- units::set_units(totshp$area, "acre")
    sf <- totshp %>% 
      st_join(sf %>%
                mutate(uniquegeom = ifelse(lengths(st_equals(sf$geometry))>1, 1, 0)) %>%
                filter(uniquegeom==0) %>%
                dplyr::select(-uniquegeom) %>%
                st_make_valid()) %>%
      rename(acresid = id)  
    
  }else {
    
    sf <- sf %>%
      mutate(AcresMade = units::set_units(st_area(geometry), "acre"))
    
  }
  
  
  ##### Download SSURGO data
  sf.areas <- get_ssurgo_fixed(
    template = sf,
    label = title, 
    force.redo = TRUE
  )
  # NOTE: get_ssurgo_fixed automatically restores original FedData functions
  # after running, so subsequent FedData calls (e.g. get_ghcn_daily) are unaffected
  
  sf <- st_transform(sf, 4326)
  sfsave <- st_transform(sfsave, 4326)
  
  ### Get Soil Temperature
  ## https://zenodo.org/records/4558732
  soiltemps.r <- rast("~/OneDrive - The Nature Conservancy/CPRG/Data/SBIO1_Annual_Mean_Temperature_5_15cm.tif")
  sf$soiltemp <- exact_extract(soiltemps.r, sf, "weighted_mean", weights="area")
  
  ### Clean spatial data to join with other tables
  sf.spatial <- sf.areas$spatial %>%
    mutate(MUKEY = as.numeric(MUKEY))
  colnames(sf.spatial) <- tolower(colnames(sf.spatial))
  
  #### Build master database 
  sf.sub.surg <- sf.spatial %>%
    dplyr::select_if(not_all_na) %>%
    mutate(numacres = as.vector(st_area(geom) * 0.000247105)) %>%
    distinct() %>%
    ### Get soil types to assess muck and water table depth
    left_join(sf.areas$tabular$muaggatt %>% dplyr::select(muname, musym, wtdepannmin, wtdepaprjunmin),
              relationship = "many-to-many") %>%
    ## Get unique ids and vegetation type
    left_join(sf.areas$tabular$component %>% dplyr::select(mukey, cokey, earthcovkind1) %>% distinct(),
              relationship = "many-to-many") %>%
    ## Get Bulk Density
    left_join(sf.areas$tabular$chorizon %>% dplyr::select(cokey, dbovendry.r) %>% 
                rename(bulkdensity = dbovendry.r)) %>%
    ## Add in industrial forest layers
    left_join(sf.areas$tabular$coforprod %>% dplyr::select(cokey, siteindex.r) %>% distinct() %>% 
                mutate(siteindex.r = ifelse(is.na(siteindex.r), 0, siteindex.r)),
              relationship = "many-to-many") %>%
    distinct() %>%
    na.omit() %>%
    dplyr::select(-cokey) %>%
    group_by(musym, mukey, geom, muname, wtdepannmin, wtdepaprjunmin, earthcovkind1, bulkdensity, numacres) %>%
    summarize(industrial = mean(siteindex.r)) %>%
    ungroup()
  
  sf.sub.surg <- st_intersection(sf.sub.surg, sf)  
  
  ### Build final GHG predictions mod
  sfindex <- sf.sub.surg %>%
    ungroup() %>%
    st_make_valid() %>%
    rename(muckname = muname,
           watertable = wtdepannmin) %>%
    mutate(mucktype = ifelse(grepl("muck", muckname)==TRUE, 1, 0)) %>%
    group_by(acresid) %>%
    summarize(vegetation = earthcovkind1[which.max(numacres)], 
              watertable = max(watertable, na.rm=TRUE),
              bulkdensity = mean(bulkdensity, na.rm=TRUE),
              mucktype = max(mucktype, na.rm=TRUE),
              soiltemp = mean(soiltemp, na.rm=TRUE),
              numacres = sum(unique(numacres), na.rm=TRUE)) %>%
    ungroup() 
  
  
  sfoutput <- list(sf.sub.surg, sfindex, sf, sfsave)
  
  return(sfoutput)
  
}