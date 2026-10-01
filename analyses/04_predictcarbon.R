#### Update Carbon Emissions Model and test results
### Started 24 October 2025 by Cat

### 14 May 2026 - adding in more sites

### housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## Load Libraries
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyterra)
library(tidyr)
library(terra)
library(patchwork)
library(broom.mixed)
library(sf)
library(brms)

### Set working directory
setwd("~/Documents/git/tnc_cprg/analyses/peatlands/")

#### Load in Shapefile
sf <- st_read("../../data/North Carolina/NC_Pocosin_Restoration_Sites_2026/")

## Load in full dataset to get z values
cleandat <- read.csv("output/clean_swails.csv")

### Load in models
load("models/soilresp.Rdata")
load("models/soilresp.simple.Rdata")
load("models/soilresp.simple.wtd.Rdata")

# Define the sites with their names, filenames, and filter criteria
sites <- list(
  list(name = "Van Swamp", filename = "vanswamp", proj_name = "Van Swamp Restoration"),
  list(name = "Hofmann Forest", filename = "hofmannforest", proj_name = "Hofmann Forest Pocosin Rewetting"),
  list(name = "Holly Shelter", filename = "hollyshelter", proj_name = "Holly Shelter Pocosin"),
  list(name = "Angola Bay", filename = "angolabay", proj_name = "Angola Bay Restoration Area")
)

# Define the years
#years <- c(2008, 2012, 2016, 2020, 2023)
years <- c(2008, 2012, 2015, 2016, 2020, 2023)

# Loop through each site
for (site in sites) {
  # Filter the shapefile for this site
  sf_site <- sf %>%
    filter(Proj_Name == site$proj_name)
  
  title <- site$name
  filename_base <- site$filename
  
  # Loop through each year
  for (inputyear in years) {
    filename <- paste0(filename_base, "_", inputyear)
    
    ### Load in restored data
    carbonestimates <- read.csv(paste0("output/clean_", filename,"_restored.csv")) %>%
      mutate(wtd = ifelse(tx == 1, 10, wtd)) ## update to within 20cm of soil surface
    
    ### Get z scores for parameters
    moddat <- cleandat %>%
      full_join(carbonestimates %>% mutate(wtd.old = wtd)) %>%
      mutate(wtd.z = (wtd - mean(wtd, na.rm=TRUE)) / (2 * sd(wtd, na.rm=TRUE)),
             wtd.old.z = (wtd.old - mean(wtd.old, na.rm=TRUE)) / (2 * sd(wtd.old, na.rm=TRUE)),
             soiltemp.z = (soiltemp - mean(soiltemp, na.rm=TRUE)) / (2 * sd(soiltemp, na.rm=TRUE)),
             pdsi.z = (scPDSI - mean(scPDSI, na.rm=TRUE)) / (2 * sd(scPDSI, na.rm=TRUE)),
             gdd.z = (totalGDD - mean(totalGDD, na.rm=TRUE)) / (2 * sd(totalGDD, na.rm=TRUE)),
             elev.z = (elevation - mean(elevation, na.rm=TRUE)) / (2 * sd(elevation, na.rm=TRUE)),
             dist.z = (distkm - mean(distkm, na.rm=TRUE)) / (2 * sd(distkm, na.rm=TRUE)),
             tx.z = (tx - mean(tx, na.rm=TRUE)) / (2 * sd(tx, na.rm=TRUE)),
             year.z = (year - mean(year, na.rm=TRUE)) / (2 * sd(year, na.rm=TRUE)),
             lai.z = (lai - mean(lai, na.rm=TRUE)) / (2 * sd(lai, na.rm=TRUE)),
             ndvi.z = (ndvi - mean(ndvi, na.rm=TRUE)) / (2 * sd(ndvi, na.rm=TRUE)),
      )
    
    testsite <- moddat %>%
      filter(location == title)
    
    testsite$soilresp_slopes <- predict(soilrespmod, newdat = testsite, allow_new_levels = TRUE)[,"Estimate"]
    testsite$soilresp.error_slopes <- predict(soilrespmod, newdat = testsite,
                                              probs = c(0.11, 0.89), allow_new_levels = TRUE)[,"Est.Error"]
    
    
    #### Build model replicating original model
    ### Ensure the same outputs every time
    set.seed(1221)
    testsite$soilresp.simp <- predict(soilresp.simple, newdat = testsite, allow_new_levels = TRUE)[,"Estimate"]
    testsite$soilresp.simp.error <- predict(soilresp.simple, newdat = testsite, 
                                            probs = c(0.11, 0.89), allow_new_levels = TRUE)[,"Est.Error"]
    
    #### Build model replicating original model w. just WTD
    testsite$soilresp.simp.wtd <- predict(soilresp.simple.wtd, newdat = testsite, allow_new_levels = TRUE)[,"Estimate"]
    testsite$soilresp.simp.wtd.error <- predict(soilresp.simple.wtd, newdat = testsite, 
                                            probs = c(0.11, 0.89), allow_new_levels = TRUE)[,"Est.Error"]
    
    
    #### Now for the original Swails model
    originalmod <- lm(soilresp ^ 0.25 ~ wtd + soiltemp, data=cleandat)
    
    originalmod.pred <- predict(originalmod, newdata = testsite, se.fit = TRUE)
    
    testsite$swails.soilresp.error <- sqrt(originalmod.pred$se.fit^2 + originalmod.pred$residual.scale^2)

    
    clean_estimates <- testsite %>%
      mutate(
        ### Final model
        totalc_slopes = (soilresp_slopes ^ 4) * 0.5 * 35.45,
        totalc_slopes_highratio = (soilresp_slopes ^ 4) * 0.7 * 35.45,
        totalc_slopes_lowratio = (soilresp_slopes ^ 4) * 0.3 * 35.45,
        totalc_slopes_mixedratio = ifelse(tx == 1, (soilresp_slopes ^ 4) * 0.5 * 35.45, (soilresp_slopes ^ 4) * 0.7 * 35.45),
        ### Simple comparison using Bayesian
        totalc.simp = (soilresp.simp ^ 4) * 0.5 * 35.45,
        totalc.simp_highratio = (soilresp.simp ^ 4) * 0.7 * 35.45,
        totalc.simp_lowratio = (soilresp.simp ^ 4) * 0.3 * 35.45,
        totalc.simp_mixedratio = ifelse(tx == 1, (soilresp.simp ^ 4) * 0.5 * 35.45, (soilresp.simp ^ 4) * 0.7 * 35.45),
        ### Just WTD comparison using Bayesian
        totalc.simp.wtd = (soilresp.simp.wtd ^ 4) * 0.5 * 35.45,
        totalc.simp.wtd_highratio = (soilresp.simp.wtd ^ 4) * 0.7 * 35.45,
        totalc.simp.wtd_lowratio = (soilresp.simp.wtd ^ 4) * 0.3 * 35.45,
        totalc.simp.wtd_mixedratio = ifelse(tx == 1, (soilresp.simp.wtd ^ 4) * 0.5 * 35.45, (soilresp.simp.wtd ^ 4) * 0.7 * 35.45),
        ### Swails original
        swails.totalc = ((coef(originalmod)[1] + coef(originalmod)[2] * wtd
                                   + coef(originalmod)[3] * soiltemp) ^ 4) * 35.45 * 0.5,
        swails.totalc_highratio = ((coef(originalmod)[1] + coef(originalmod)[2] * wtd + 
                                               coef(originalmod)[3] * soiltemp) ^ 4) * 35.45 * 0.7,
        swails.totalc_lowratio = ((coef(originalmod)[1] + coef(originalmod)[2] * wtd +
                                              coef(originalmod)[3] * soiltemp) ^ 4) * 35.45 * 0.3,
        swails.totalc_mixedratio = ifelse(tx == 1, ((coef(originalmod)[1] + coef(originalmod)[2] * wtd +
                                                                coef(originalmod)[3] * soiltemp) ^ 4) 
                                                   * 35.45 * 0.5,((coef(originalmod)[1] + coef(originalmod)[2] * wtd +
                                                                     coef(originalmod)[3] * soiltemp) ^ 4) 
                                                   * 35.45 * 0.7),
      )
    
    estimates_wide <- clean_estimates %>%
      st_drop_geometry() %>%
      select(site, tx, totalc_slopes:swails.totalc_mixedratio, acresid, numacres) %>%
      mutate(tx = ifelse(tx == 0, "drained", "restored")) %>%
      pivot_wider(
        names_from = tx,
        values_from = c(totalc_slopes, totalc_slopes_highratio, totalc_slopes_lowratio, totalc_slopes_mixedratio, 
                        totalc.simp, totalc.simp_highratio, totalc.simp_lowratio, totalc.simp_mixedratio, 
                        totalc.simp.wtd, totalc.simp.wtd_highratio, totalc.simp.wtd_lowratio, totalc.simp.wtd_mixedratio, 
                        swails.totalc, swails.totalc_highratio, swails.totalc_lowratio, swails.totalc_mixedratio,)
      ) %>%
      mutate(### Final Model
        addit_slopes = (totalc_slopes_drained - totalc_slopes_restored),
        addit_slopes_highratio = (totalc_slopes_highratio_drained - totalc_slopes_highratio_restored),
        addit_slopes_lowratio = (totalc_slopes_lowratio_drained - totalc_slopes_lowratio_restored),
        addit_slopes_mixedratio = (totalc_slopes_mixedratio_drained - totalc_slopes_mixedratio_restored),
        ### Simple comparison
        addit.simp = (totalc.simp_drained - totalc.simp_restored),
        addit.simp_highratio = (totalc.simp_highratio_drained - totalc.simp_highratio_restored),
        addit.simp_lowratio = (totalc.simp_lowratio_drained - totalc.simp_lowratio_restored),
        addit.simp_mixedratio = (totalc.simp_mixedratio_drained - totalc.simp_mixedratio_restored),
        ### Just WTD comparison
        addit.simp.wtd = (totalc.simp.wtd_drained - totalc.simp.wtd_restored),
        addit.simp.wtd_highratio = (totalc.simp.wtd_highratio_drained - totalc.simp.wtd_highratio_restored),
        addit.simp.wtd_lowratio = (totalc.simp.wtd_lowratio_drained - totalc.simp.wtd_lowratio_restored),
        addit.simp.wtd_mixedratio = (totalc.simp.wtd_mixedratio_drained - totalc.simp.wtd_mixedratio_restored),
        ### Original Swails 
        addit.swails = (swails.totalc_drained - swails.totalc_restored),
        addit.swails_highratio = (swails.totalc_highratio_drained - swails.totalc_highratio_restored),
        addit.swails_lowratio = (swails.totalc_lowratio_drained - swails.totalc_lowratio_restored),
        addit.swails_mixedratio = (swails.totalc_mixedratio_drained - swails.totalc_mixedratio_restored))
    
    all_error <- clean_estimates %>%
      st_drop_geometry() %>%
      select(site, tx, acresid, numacres, soilresp.error_slopes, 
             soilresp.simp.error, soilresp.simp.wtd.error,
             swails.soilresp.error) %>%
      mutate(tx = ifelse(tx == 0, "drained", "restored"))
    
    estimates_wide <- left_join(estimates_wide, all_error)
    
    totalacres <- sf_site %>%
      mutate(acres = as.vector(st_area(geometry) * 0.000247105))
    
    results = estimates_wide %>%
      summarize(
        # Final Model
        avg_addit_slopes = weighted.mean(addit_slopes, numacres),
        err_addit_slopes = weighted.mean(soilresp.error_slopes, numacres),
        avg_addit_slopes_highratio = weighted.mean(addit_slopes_highratio, numacres),
        err_addit_slopes_highratio = weighted.mean(soilresp.error_slopes, numacres),
        avg_addit_slopes_lowratio = weighted.mean(addit_slopes_lowratio, numacres),
        err_addit_slopes_lowratio = weighted.mean(soilresp.error_slopes, numacres),
        avg_addit_slopes_mixedratio = weighted.mean(addit_slopes_mixedratio, numacres),
        err_addit_slopes_mixedratio = weighted.mean(soilresp.error_slopes, numacres),
        # Simple comparison
        avg_additsimp = weighted.mean(addit.simp, numacres),
        err_additsimp = weighted.mean(soilresp.simp.error, numacres),
        avg_additsimp_highratio = weighted.mean(addit.simp_highratio, numacres),
        err_additsimp_highratio = weighted.mean(soilresp.simp.error, numacres),
        avg_additsimp_lowratio = weighted.mean(addit.simp_lowratio, numacres),
        err_additsimp_lowratio = weighted.mean(soilresp.simp.error, numacres),
        avg_additsimp_mixedratio = weighted.mean(addit.simp_mixedratio, numacres),
        err_additsimp_mixedratio = weighted.mean(soilresp.simp.error, numacres),
        # Just WTD comparison
        avg_additsimp.wtd = weighted.mean(addit.simp.wtd, numacres),
        err_additsimp.wtd = weighted.mean(soilresp.simp.wtd.error, numacres),
        avg_additsimp.wtd_highratio = weighted.mean(addit.simp.wtd_highratio, numacres),
        err_additsimp.wtd_highratio = weighted.mean(soilresp.simp.wtd.error, numacres),
        avg_additsimp.wtd_lowratio = weighted.mean(addit.simp.wtd_lowratio, numacres),
        err_additsimp.wtd_lowratio = weighted.mean(soilresp.simp.wtd.error, numacres),
        avg_additsimp.wtd_mixedratio = weighted.mean(addit.simp.wtd_mixedratio, numacres),
        err_additsimp.wtd_mixedratio = weighted.mean(soilresp.simp.wtd.error, numacres),
        # Original Swails
        swails.avg_addit = weighted.mean(addit.swails, numacres),
        swails.err_addit = weighted.mean(swails.soilresp.error, numacres),
        swails.avg_addit_highratio = weighted.mean(addit.swails_highratio, numacres),
        swails.err_addit_highratio = weighted.mean(swails.soilresp.error, numacres),
        swails.avg_addit_lowratio = weighted.mean(addit.swails_lowratio, numacres),
        swails.err_addit_lowratio = weighted.mean(swails.soilresp.error, numacres),
        swails.avg_addit_mixedratio = weighted.mean(addit.swails_mixedratio, numacres),
        swails.err_addit_mixedratio = weighted.mean(swails.soilresp.error, numacres)
      )
    
    ################################################################################
    ################################################################################
    ### Save the summary results to a dataframe and then plot figure and table
    
    ### Compute shared limits from the data
    shared_ylim <- range(c(results$avg_addit_slopes, results$avg_additsimp,
                           results$avg_addit_slopes_highratio, results$avg_additsimp_highratio,
                           results$avg_addit_slopes_lowratio, results$avg_additsimp_lowratio,
                           results$avg_addit_slopes_mixedratio, results$avg_additsimp_mixedratio,
                           results$avg_additsimp.wtd, results$avg_additsimp.wtd_highratio,
                           results$avg_additsimp.wtd_lowratio, results$avg_additsimp.wtd_mixedratio), 
                         na.rm = TRUE)
    
    ### Pivot to long format for global average
    results_long <- data.frame(
      model    = c("Updated model",
                   "Original model", "Original model - WTD only"),
      estimate = c(results$avg_addit_slopes, results$avg_additsimp, results$avg_additsimp.wtd),
      se       = c(results$err_addit_slopes, results$err_additsimp, results$err_additsimp.wtd)
    ) %>%
      mutate(lower = estimate - se,
             upper = estimate + se,
             model = factor(model, levels = c("Original model - WTD only",
                                              "Original model",
                                              "Updated model")),
             ratio = "original")
    
    ### Pivot to long format for high ratio
    results_long_highratio <- data.frame(
      model    = c("Updated model",
                   "Original model", "Original model - WTD only"),
      estimate = c(results$avg_addit_slopes_highratio, results$avg_additsimp_highratio, results$avg_additsimp.wtd_highratio),
      se       = c(results$err_addit_slopes_highratio, results$err_additsimp_highratio, results$err_additsimp.wtd_highratio)
    ) %>%
      mutate(lower = estimate - se,
             upper = estimate + se,
             model = factor(model, levels = c("Original model - WTD only", "Original model",
                                              "Updated model")),
             ratio = "high")
    
    ### Pivot to long format for low ratio
    results_long_lowratio <- data.frame(
      model    = c("Updated model",
                   "Original model", "Original model - WTD only"),
      estimate = c(results$avg_addit_slopes_lowratio, results$avg_additsimp_lowratio, results$avg_additsimp.wtd_lowratio),
      se       = c(results$err_addit_slopes_lowratio, results$err_additsimp_lowratio, results$err_additsimp.wtd_lowratio)
    ) %>%
      mutate(lower = estimate - se,
             upper = estimate + se,
             model = factor(model, levels = c("Original model - WTD only", "Original model",
                                              "Updated model")),
             ratio = "low")
    
    ### Pivot to long format for mixed ratio
    results_long_mixedratio <- data.frame(
      model    = c("Updated model",
                   "Original model", "Original model - WTD only"),
      estimate = c(results$avg_addit_slopes_mixedratio, results$avg_additsimp_mixedratio, results$avg_additsimp.wtd_mixedratio),
      se       = c(results$err_addit_slopes_mixedratio, results$err_additsimp_mixedratio, results$err_additsimp.wtd_mixedratio)
    ) %>%
      mutate(lower = estimate - se,
             upper = estimate + se,
             model = factor(model, levels = c("Original model - WTD only", "Original model",
                                              "Updated model"
             )),
             ratio = "mixed")
    
    results_long <- full_join(results_long, results_long_highratio) %>%
      full_join(results_long_lowratio) %>%
      full_join(results_long_mixedratio)
    
    ### Plot
    cgains <- ggplot(results_long, 
                     aes(x = model, y = estimate, 
                         color = ratio,
                         group = interaction(model, ratio))) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
      geom_pointrange(aes(ymin = lower, ymax = upper),
                      size     = 0.4, 
                      linewidth = 1) +
      scale_color_manual(
        values = c("original" = "#332288", 
                   "high" = "#E66101", 
                   "low" = "#1B9E77", 
                   "mixed" = "#B35806"),
        labels = c("original" = "50% (Global Average)",
                   "high"  = "70% (High end of range)",
                   "low"  = "30% (Low end of range)",
                   "mixed"  = "50% for restored, 70% for drained")) +
      scale_x_discrete(labels = c(
        "Original model" = "Original model",
        "Updated model" = "Updated model",
        "Original model - WTD only" = "Original model - WTD only"
      )) +
      coord_cartesian(ylim = c(0, 7)) +
      labs(x = NULL, y = expression("t CO"[2]*"e acre"^{-1}*" yr"^{-1}),
        title = paste0(title, " ", inputyear, " - Carbon benefit estimates"),
        subtitle = "Point = weighted mean; bars = 89% CI",
        color = "Heterotrophic\npartitioning ratio"
      ) +
      theme_bw() +
      theme(legend.position = "right", plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9, color = "gray40"), 
        axis.text.x = element_text(size = 9), panel.grid.major.x = element_blank())
    
    
    grid.arrange(results_long %>% 
                   arrange(model) %>%
                   select(-upper, -lower) %>%
                   mutate(ratio = ifelse(ratio == "original", "50% (Global average)", 
                                         ifelse(ratio == "highratio", "70% (High end of range)",
                                                ifelse(ratio == "lowratio", 
                                                       "30% (Low end of range", "50% for restored, 70% for drained"))),
                          estimate = round(estimate, digits = 2),
                          se = round(se, digits = 2)) %>%
                   rename(Model = model,
                          AnnualCarbon = estimate, 
                          Error = se,
                          Ratio = ratio) %>%
                   tableGrob(theme = ttheme_default(
                     core = list(bg_params=list(fill=c("grey90", "white"))),
                     colhead = list(fg_params=list(col="white"),
                                    bg_params=list(fill="green4"))), rows = NULL))
    
    
    
    # Save figure
    png(paste0("figures/fullmodelcompare_figure_", filename, ".png"), 
        width=7,
        height=6, units="in", res = 350 )
    print(cgains)
    dev.off()
    
    # Save table
    png(paste0("figures/fullmodelcompare_table_", filename, ".png"), 
        width=7,
        height=6, units="in", res = 350 )
    grid.arrange(results_long %>% 
                   arrange(model) %>%
                   select(-upper, -lower) %>%
                   mutate(ratio = ifelse(ratio == "original", "50% (Global average)", 
                                         ifelse(ratio == "highratio", "70% (High end of range)",
                                                ifelse(ratio == "lowratio", 
                                                       "30% (Low end of range", "50% for restored, 70% for drained"))),
                          estimate = round(estimate, digits = 2),
                          se = round(se, digits = 2)) %>%
                   rename(Model = model,
                          AnnualCarbon = estimate, 
                          Error = se,
                          Ratio = ratio) %>%
                   tableGrob(theme = ttheme_default(
                     core = list(bg_params=list(fill=c("grey90", "white"))),
                     colhead = list(fg_params=list(col="white"),
                                    bg_params=list(fill="green4"))), rows = NULL))
    dev.off()
    
    write.csv(results_long, paste0("output/carbonestimates_", filename,".csv"), row.names = FALSE)
    
    print(paste0("Completed ", title, " for year ", inputyear))
    
  } # end year loop
  
} # end site loop
