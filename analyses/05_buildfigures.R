### Plot all outputs into consolidated figures
### Started 15 May 2026 by Cat

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
library(sf)
library(forcats)
library(flextable)
library(webshot2)
library(brms)
library(bayesplot)
library(broom.mixed)

### Ensure the same outputs every time
set.seed(1221)


### Set working directory
setwd("~/Documents/git/tnc_cprg/analyses/peatlands/")

### Load in models
#load("models/totalc.Rdata")
load("models/soilresp.Rdata")
load("models/soilresp.noelev.Rdata")
load("models/soilresp_study.Rdata")

## Load in full dataset to get z values
cleandat <- read.csv("output/clean_swails.csv")


### Load all dataframes
van <- read.csv("output/carbonestimates_vanswamp_2008.csv") %>%
  mutate(year = 2008, name = "Van Swamp") %>%
  rbind(read.csv("output/carbonestimates_vanswamp_2012.csv") %>%
          mutate(year = 2012, name = "Van Swamp")) %>%
  rbind(read.csv("output/carbonestimates_vanswamp_2016.csv") %>%
          mutate(year = 2016, name = "Van Swamp")) %>%
  rbind(read.csv("output/carbonestimates_vanswamp_2020.csv") %>%
          mutate(year = 2020, name = "Van Swamp")) %>%
  rbind(read.csv("output/carbonestimates_vanswamp_2023.csv") %>%
          mutate(year = 2023, name = "Van Swamp"))

hofmann <- read.csv("output/carbonestimates_hofmannforest_2008.csv") %>%
  mutate(year = 2008, name = "Hofmann Forest") %>%
  rbind(read.csv("output/carbonestimates_hofmannforest_2012.csv") %>%
          mutate(year = 2012, name = "Hofmann Forest")) %>%
  rbind(read.csv("output/carbonestimates_hofmannforest_2016.csv") %>%
          mutate(year = 2016, name = "Hofmann Forest")) %>%
  rbind(read.csv("output/carbonestimates_hofmannforest_2020.csv") %>%
          mutate(year = 2020, name = "Hofmann Forest")) %>%
  rbind(read.csv("output/carbonestimates_hofmannforest_2023.csv") %>%
          mutate(year = 2023, name = "Hofmann Forest"))

holly <- read.csv("output/carbonestimates_hollyshelter_2008.csv") %>%
  mutate(year = 2008, name = "Holly Shelter") %>%
  rbind(read.csv("output/carbonestimates_hollyshelter_2015.csv") %>%
          mutate(year = 2015, name = "Holly Shelter")) %>%
  rbind(read.csv("output/carbonestimates_hollyshelter_2016.csv") %>%
          mutate(year = 2016, name = "Holly Shelter")) %>%
  rbind(read.csv("output/carbonestimates_hollyshelter_2020.csv") %>%
          mutate(year = 2020, name = "Holly Shelter")) %>%
  rbind(read.csv("output/carbonestimates_hollyshelter_2023.csv") %>%
          mutate(year = 2023, name = "Holly Shelter"))


angola <- read.csv("output/carbonestimates_angolabay_2008.csv") %>%
  mutate(year = 2008, name = "Angola Bay") %>%
  rbind(read.csv("output/carbonestimates_angolabay_2012.csv") %>%
          mutate(year = 2012, name = "Angola Bay")) %>%
  rbind(read.csv("output/carbonestimates_angolabay_2016.csv") %>%
          mutate(year = 2016, name = "Angola Bay")) %>%
  rbind(read.csv("output/carbonestimates_angolabay_2020.csv") %>%
          mutate(year = 2020, name = "Angola Bay")) %>%
  rbind(read.csv("output/carbonestimates_angolabay_2023.csv") %>%
          mutate(year = 2023, name = "Angola Bay"))



##### Bring them all together
all <- bind_rows(list(van, hofmann, holly, angola)) %>%
  filter(!model %in% c("Simple WTD\n(Bayesian - Swails Comparison)",
                       "Original model - WTD only")) %>%
  mutate(name = case_when(
           name == "Van Swamp" ~ "Site VS",
           name == "Hofmann Forest" ~ "Site HF",
           name == "Holly Shelter" ~ "Site HS",
           name == "Angola Bay" ~ "Site AB"))

annual.p <- ggplot(all, aes(x = year, y = estimate, col = ratio)) + #, linetype = model
  geom_point() + geom_line() + theme_bw() +
  scale_color_manual(
    values = c("original" = "#332288", 
               "high" = "#E66101", 
               "low" = "#1B9E77", 
               "mixed" = "#B35806"),
    labels = c("original" = "50% (Global Average)",
               "high"  = "70% (High end of range)",
               "low"  = "30% (Low end of range)",
               "mixed"  = "50% for restored, 70% for drained")) +
  #scale_linetype_manual(
  #  values = c("Original model" = "solid",
  #             "Updated model" = "dashed")) +
  facet_grid(name~model) +
  labs(x = NULL, y = expression("Mg CO"[2]*" acre"^{-1}*" yr"^{-1}),
    title = "Annual carbon benefit estimates",
    color = "Heterotrophic\npartitioning ratio", 
    linetype = "Model") +
  theme_bw() +
  theme(legend.position = "right", plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray40"), 
    axis.text.x = element_text(size = 9), panel.grid.major.x = element_blank())

##### Try condensing the figure rather than breaking out by year
site.summary <- all %>%
  group_by(name, model, ratio) %>%
  summarise(
    mean = mean(estimate, na.rm = TRUE),
    se = sd(estimate, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  mutate(ratio = ifelse(ratio == "original", "global", ratio),
         model = ifelse(model == "Original model", "Baseline model", model))

### Save as table
summary_table <- site.summary %>%
  pivot_wider(
    names_from = ratio,
    values_from = c(mean, se)
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 2))) %>%
  arrange(name, model)

# Add in breaks
site_breaks <- which(summary_table$name[-1] != summary_table$name[-nrow(summary_table)])

## Double check it is working
site_breaks

summary_table %>%
  flextable() %>%
  set_header_labels(name = "Site", model = "Model") %>%
  theme_booktabs() %>%
  hline(i = site_breaks, part = "body", border = officer::fp_border(color = "gray80",width = 1)) %>%
  autofit() %>%
  save_as_image(path = "figures/table_avggains_compare.png", width = 7, height = 6)

### Find differences for MS
site.summary %>%
  pivot_wider(names_from = ratio,
              values_from = c(mean, se)) %>%
  mutate(across(where(is.numeric), ~round(., 2))) %>%
  mutate(diff = mean_mixed - mean_low) %>%
  #group_by(model) %>%
  summarize(avgdiff = mean(diff),
            sddiff = sd(diff))

site.summary %>%
  pivot_wider(names_from = model,
              values_from = c(mean, se)) %>%
  mutate(across(where(is.numeric), ~round(., 2))) %>%
  mutate(diff = `mean_Updated model` - `mean_Baseline model`) %>%
  #group_by(model) %>%
  summarize(avgdiff = mean(diff),
            sddiff = sd(diff))


site.p <- ggplot(site.summary, 
                 aes(x = name, y = mean, 
                     color = ratio, 
                     group = interaction(ratio, model))) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, linewidth = 0.8) +
  scale_color_manual(
    values = c("original" = "#332288", 
               "high" = "#E66101", 
               "low" = "#1B9E77", 
               "mixed" = "#B35806"),
    labels = c("original" = "50% (Global Average)",
               "high"  = "70% (High end of range)",
               "low"  = "30% (Low end of range)",
               "mixed"  = "50% for restored, 70% for drained")) +
  facet_wrap(~model) +
  labs(x = "", y = expression("Mg CO"[2]*" acre"^{-1}*" yr"^{-1}),
    title = "Mean annual carbon benefit estimates by site",
    subtitle = "Error bars = standard error across years",
    color = "Heterotrophic\npartitioning ratio") +
  theme_bw() +
  theme(legend.position = "right", plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    panel.grid.major.x = element_blank())

png(paste0("figures/avggains_compare.png"), 
    width=7,
    height=4, units="in", res = 350 )
site.p
dev.off()

##### Look at total over time
sf <- st_read("../../data/North Carolina/NC_Pocosin_Restoration_Sites_2026/") %>%
  mutate(name = case_when(
    Proj_Name == "Van Swamp Restoration" ~ "Site VS",
    Proj_Name == "Hofmann Forest Pocosin Rewetting" ~ "Site HF",
    Proj_Name == "Holly Shelter Pocosin" ~ "Site HS",
    Proj_Name == "Angola Bay Restoration Area" ~ "Site AB"),
    acres = as.vector(st_area(geometry) * 0.000247105))

grouped <- left_join(all, sf %>% select(name, acres) %>% st_drop_geometry) %>%
  mutate(model = ifelse(model == "Original model", "Baseline model", model)) %>%
  group_by(model, ratio, name) %>%
  arrange(year) %>%
  mutate(timediff = year - lag(year, default = first(year))) %>%
  ungroup() %>%
  mutate(totalc = estimate * acres * timediff) %>%
  group_by(name, model, ratio) %>%
  mutate(totalc = ifelse(is.na(totalc), totalc[year==2008], totalc),
         totalc = cumsum(totalc),
         ratio = ifelse(ratio == "original", "50% (Global average)", 
                        ifelse(ratio == "high", "70% (High end of range)",
                               ifelse(ratio == "low", 
                                      "30% (Low end of range)", 
                                      "50% for restored, 70% for drained"))))

# Reorder ratio by maximum totalc from Updated model
max_by_ratio <- grouped %>%
  filter(model == "Updated model", name == "Site HF", year == 2023) %>%
  group_by(ratio) %>%
  summarise(max_totalc = max(totalc, na.rm = TRUE), .groups = 'drop')

grouped <- grouped %>%
  left_join(max_by_ratio, by = "ratio") %>%
  mutate(ratio = fct_reorder(ratio, max_totalc))




total.p <- ggplot(grouped, aes(x = year, y = totalc/1000000, col = ratio, linetype = model)) +
  geom_point() + 
  geom_line() + 
  scale_color_manual(
    values = c("50% (Global average)" = "#332288", 
               "70% (High end of range)" = "#E66101", 
               "30% (Low end of range)" = "#1B9E77", 
               "50% for restored, 70% for drained" = "#B35806")) +
  scale_linetype_manual(
    values = c("Baseline model" = "solid",
               "Updated model" = "dashed",
               "Baseline model (WTD only)" = "dotted")) +
  facet_grid(name~ratio) +
  labs(x = NULL, y = expression("MMg CO"[2]),
    title = "Accumulated carbon benefit estimates",
    color = "Heterotrophic\npartitioning ratio",
    linetype = "Model") +
  theme_bw() +
  theme(legend.position = "right", plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    axis.text.x = element_text(size = 9), panel.grid.major.x = element_blank())

png(paste0("figures/accumulatedgains_compare.png"), 
    width=12,
    height=6, units="in", res = 350 )
total.p
dev.off()


#################################################################################
#### Show spatiotemporal variation across sites and influence on annual estimates
ab <- read.csv("output/clean_angolabay_2008_restored.csv") %>% mutate(year = 2008) %>%
  full_join(read.csv("output/clean_angolabay_2012_restored.csv") %>% mutate(year = 2012)) %>%
  full_join(read.csv("output/clean_angolabay_2016_restored.csv") %>% mutate(year = 2016)) %>%
  full_join(read.csv("output/clean_angolabay_2020_restored.csv") %>% mutate(year = 2020)) %>%
  full_join(read.csv("output/clean_angolabay_2023_restored.csv") %>% mutate(year = 2023))

vs <- read.csv("output/clean_vanswamp_2008_restored.csv") %>% mutate(year = 2008) %>%
  full_join(read.csv("output/clean_vanswamp_2012_restored.csv") %>% mutate(year = 2012)) %>%
  full_join(read.csv("output/clean_vanswamp_2016_restored.csv") %>% mutate(year = 2016)) %>%
  full_join(read.csv("output/clean_vanswamp_2020_restored.csv") %>% mutate(year = 2020)) %>%
  full_join(read.csv("output/clean_vanswamp_2023_restored.csv") %>% mutate(year = 2023))

hf <- read.csv("output/clean_hofmannforest_2008_restored.csv") %>% mutate(year = 2008) %>%
  full_join(read.csv("output/clean_hofmannforest_2012_restored.csv") %>% mutate(year = 2012)) %>%
  full_join(read.csv("output/clean_hofmannforest_2016_restored.csv") %>% mutate(year = 2016)) %>%
  full_join(read.csv("output/clean_hofmannforest_2020_restored.csv") %>% mutate(year = 2020)) %>%
  full_join(read.csv("output/clean_hofmannforest_2023_restored.csv") %>% mutate(year = 2023))

hs <- read.csv("output/clean_hollyshelter_2008_restored.csv") %>% mutate(year = 2008) %>%
  full_join(read.csv("output/clean_hollyshelter_2015_restored.csv") %>% mutate(year = 2015)) %>%
  full_join(read.csv("output/clean_hollyshelter_2016_restored.csv") %>% mutate(year = 2016)) %>%
  full_join(read.csv("output/clean_hollyshelter_2020_restored.csv") %>% mutate(year = 2020)) %>%
  full_join(read.csv("output/clean_hollyshelter_2023_restored.csv") %>% mutate(year = 2023))


##### Bring them all together
all_preds <- bind_rows(list(vs, hf, hs, ab)) %>%
  mutate(site = case_when(
    site == "Van Swamp" ~ "Site VS",
    site == "Hofmann Forest" ~ "Site HF",
    site == "Holly Shelter" ~ "Site HS",
    site == "Angola Bay" ~ "Site AB")
  )

#### Visualize predictor variability across sites and years
climatic_variability <- all_preds %>%
  select(site, year, totalGDD, scPDSI) %>%
  pivot_longer(
    cols = c(totalGDD, scPDSI),
    names_to = "predictor",
    values_to = "value"
  ) %>%
  mutate(predictor = case_when(
    predictor == "totalGDD" ~ "Growing degree days",
    predictor == "scPDSI" ~ "Drought index (PDSI)"
  )) %>%
  distinct()

geographic_variability <- all_preds %>%
  select(site, distkm, elevation) %>%
  pivot_longer(
    cols = c(distkm, elevation),
    names_to = "predictor",
    values_to = "value"
  ) %>%
  mutate(predictor = case_when(
    predictor == "distkm" ~ "Distance to coast",
    predictor == "elevation" ~ "Elevation"
  ))

vegetative_variability <- all_preds %>%
  select(site, year, lai, ndvi) %>%
  pivot_longer(
    cols = c(lai, ndvi),
    names_to = "predictor",
    values_to = "value"
  ) %>%
  mutate(predictor = case_when(
    predictor == "lai" ~ "LAI",
    predictor == "ndvi" ~ "NDVI"
  ))

clim.p <- ggplot(climatic_variability %>% group_by(site, year, predictor) %>% 
                   summarize(sd = sd(value),
                     value = mean(value)), 
                     aes(x = factor(year), y = value, col = site, group = site)) +
  geom_point() + geom_pointrange(aes(ymax = value + sd, ymin = value - sd)) + 
  facet_wrap(~predictor, scales = "free_y", ncol = 4) +
  scale_color_manual(
    values = c("Site VS" = "#332288",
               "Site HF" = "#E66101",
               "Site HS" = "#1B9E77",
               "Site AB" = "#B35806")) +
  labs(title = "Climatic variability across sites and years", x = "", y = "Value",
    color = "Site") +
  theme_bw() +
  theme(plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 8),
    legend.position = "right")


veg.p <-  ggplot(vegetative_variability %>% group_by(site, year, predictor) %>% 
  summarize(sd = sd(value),
            value = mean(value)), 
aes(x = factor(year), y = value, col = site, group = site)) +
  geom_point() + geom_pointrange(aes(ymax = value + sd, ymin = value - sd)) + 
  facet_wrap(~predictor, scales = "free_y", ncol = 4) +
  scale_color_manual(
    values = c("Site VS" = "#332288",
               "Site HF" = "#E66101",
               "Site HS" = "#1B9E77",
               "Site AB" = "#B35806")) +
  labs(title = "Vegetative variability across sites", x = "", y = "Value",
    color = "Site") +
  theme_bw() +
  theme(plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 8),
    legend.position = "right")

geo.p <- ggplot(geographic_variability %>% group_by(site, predictor) %>% 
                  summarize(sd = sd(value),
                            value = mean(value)), 
                aes(x = site, y = value, col = site, group = site)) +
  geom_point() + geom_pointrange(aes(ymax = value + sd, ymin = value - sd)) +   
  facet_wrap(~predictor, scales = "free_y", ncol = 4) +
  scale_color_manual(
    values = c("Site VS" = "#332288",
               "Site HF" = "#E66101",
               "Site HS" = "#1B9E77",
               "Site AB" = "#B35806")) +
  labs(title = "Geographic variability across sites", x = "", y = "Value", color = "Site") +
  theme_bw() +
  theme(plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 8),
    legend.position = "right")


# Remove individual legends and add one shared legend
clim.p_clean <- clim.p + theme(legend.position = "none")
veg.p_clean <- veg.p + theme(legend.position = "none")
geo.p_clean <- geo.p + theme(legend.position = "none")

png(paste0("figures/covariates.png"), 
    width = 10, height = 10, units = "in", res = 350)
(clim.p_clean / veg.p_clean / geo.p_clean) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")
dev.off()


#################################################################################
## Plot model output to see what variables influence soil respiration
slopes.p <- tidy(soilrespmod, effects = "fixed", conf.int = TRUE, conf.level = 0.89) %>%
  filter(term != "(Intercept)") %>%
  mutate(term = case_when(
    term == "wtd.z" ~ "Water table depth",
    term == "soiltemp.z" ~ "Soil temperature",
    term == "gdd.z" ~ "Growing degree days",
    term == "lai.z" ~ "Leaf Area Index",
    term == "elev.z" ~ "Elevation",
    term == "pdsi.z" ~ "Drought index (PDSI)",
    term == "dist.z" ~ "Distance to coast",
    term == "ndvi.z" ~ "NDVI"
  )) %>%
  ggplot(aes(x = estimate, y = reorder(term, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
                  color = "purple3", 
                  size = 0.5, linewidth = 0.9) +
  geom_text(aes(label = round(estimate, 2)),
            nudge_y = 0.3, size = 3, color = "gray30") +
  labs(title = "Predictors of soil respiration - varying slopes",
    subtitle = "Posterior median \u00b1 89% CI | All predictors standardized",
    x = expression("Effect on soil respiration"^0.25),
    y = NULL) +
  theme_bw() +
  theme(plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 8,  color = "gray40"),
    axis.text.y = element_text(size = 9),
    panel.grid.major.y = element_blank())


png(paste0("figures/modeloutputs_predictors.png"), 
    width=8,
    height=6, units="in", res = 350 )
slopes.p
dev.off()


### Get z scores for parameters
moddat <- cleandat %>%
  mutate(wtd.z = (wtd - mean(wtd, na.rm=TRUE)) / (2 * sd(wtd, na.rm=TRUE)),
         #soiltemp = ifelse(is.na(soiltemp), soiltemp.clean, soiltemp),
         soiltemp.z = (soiltemp - mean(soiltemp, na.rm=TRUE)) / (2 * sd(soiltemp, na.rm=TRUE)),
         #precip.z = (precip - mean(precip, na.rm=TRUE)) / (2 * sd(precip, na.rm=TRUE)),
         pdsi.z = (scPDSI - mean(scPDSI, na.rm=TRUE)) / (2 * sd(scPDSI, na.rm=TRUE)),
         #maxt.z = (maxt - mean(maxt, na.rm=TRUE)) / (2 * sd(maxt, na.rm=TRUE)),
         gdd.z = (totalGDD - mean(totalGDD, na.rm=TRUE)) / (2 * sd(totalGDD, na.rm=TRUE)),
         #mint.z = (mint - mean(mint, na.rm=TRUE)) / (2 * sd(mint, na.rm=TRUE)),
         elev.z = (elevation - mean(elevation, na.rm=TRUE)) / (2 * sd(elevation, na.rm=TRUE)),
         dist.z = (distkm - mean(distkm, na.rm=TRUE)) / (2 * sd(distkm, na.rm=TRUE)),
         tx.z = (tx - mean(tx, na.rm=TRUE)) / (2 * sd(tx, na.rm=TRUE)),
         year.z = (year - mean(year, na.rm=TRUE)) / (2 * sd(year, na.rm=TRUE)),
         #veg.z = (vegetation - mean(vegetation, na.rm=TRUE)) / (2 * sd(vegetation, na.rm=TRUE))
         #ndvi.z = (ndvi - mean(ndvi, na.rm=TRUE)) / (2 * sd(ndvi, na.rm=TRUE)),
         lai.z = (lai - mean(lai, na.rm=TRUE)) / (2 * sd(lai, na.rm=TRUE)),
  )


#### Build model replicating Swails to see if they would predict similar outputs
### Remove NAs
moddat2 <- cleandat %>%
  mutate(date = NULL,
         uniqueid = paste0(location, plot)) 
moddat2 <- moddat2[complete.cases(moddat2),]

#### Compare uncertainty between two models:
#soilresp.simp.mod <- brm(soilresp ^ 0.25 ~ wtd, data = moddat2)
load("models/soilresp.simple.Rdata")

# Extract and prepare sigma data for plotting
sigma_data <- bind_rows(
  as_draws_df(soilresp.simple) %>%
    select(sigma) %>%
    slice_sample(n = 1000) %>%
    mutate(model = "Baseline"),
  as_draws_df(soilrespmod) %>%
    select(sigma) %>%
    slice_sample(n = 1000) %>%
    mutate(model = "Updated")
) %>%
  group_by(model) %>%
  summarise(
    estimate = median(sigma),
    conf.low = quantile(sigma, 0.055),  # (1 - 0.89) / 2
    conf.high = quantile(sigma, 0.945),
    .groups = "drop"
  )

sigma.p <- sigma_data %>%
  ggplot(aes(x = estimate, y = reorder(model, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
                  color = "purple3", 
                  size = 0.5, linewidth = 0.9) +
  geom_text(aes(label = round(estimate, 3)),
            nudge_y = 0.25, size = 3, color = "gray30") +
  labs(title = "b) Model Comparison - Sigma",
    x = "Sigma (residual SD)",
    y = NULL) +
  theme_bw() +
  theme(plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 8, color = "gray40"),
    axis.text.y = element_text(size = 9),
    panel.grid.major.y = element_blank())


### Save combined figure
combined.p <- slopes.p / sigma.p +
  plot_layout(heights = c(3, 1)) +
  plot_annotation(
    title = "Soil Respiration Model Results",
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  )

png(paste0("figures/modeloutputs_combined.png"), 
    width = 8, height = 9, units = "in", res = 350)
combined.p
dev.off()

#### Build a map of all the sites evaluated
sf <- st_read("../../data/North Carolina/NC_Pocosin_Restoration_Sites_2026/") %>%
  filter(Proj_Name %in% c("Van Swamp Restoration", "Hofmann Forest Pocosin Rewetting",
                          "Holly Shelter Pocosin", "Angola Bay Restoration Area")) %>%
  mutate(tag = c("Site AB", "Site HF", "Site VS", "Site HS"))


statenc = tigris::states() %>%
  st_transform(crs(sf)) %>%
  filter(NAME == "North Carolina")

sites_bbox <- st_bbox(sf)
statenc_cropped <- st_crop(statenc, 
                           st_buffer(st_as_sfc(sites_bbox), 50000))  # 50km buffer

countiesnc = tigris::counties(state = "NC") %>%
  st_transform(crs(sf)) 

counties_cropped = countiesnc %>%
  st_intersection(statenc_cropped)

main_map <- ggplot() + 
  geom_sf(data = counties_cropped, fill = "grey95", linewidth = 0.6) +
  geom_sf(data = sf, aes(geometry = st_centroid(geometry), color = tag), 
          size = 4, stroke = 0.5) + 
  scale_color_manual(
    values = c("Site VS" = "#332288",
               "Site HF" = "#E66101",
               "Site HS" = "#1B9E77",
               "Site AB" = "#B35806"), 
    name = "") +
  labs(title = "Test site locations in North Carolina, USA") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    legend.position = "right",
    legend.text = element_text(size = 10),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")) + coord_sf(expand = FALSE)
  

inset_map <- ggplot() +
  geom_sf(data = countiesnc, fill = "white", color = "grey92", linewidth = 0.3) +
  geom_sf(data = sf, aes(color = tag, geometry = st_centroid(geometry))) +
  scale_color_manual(
    values = c("Site VS" = "#332288",
               "Site HF" = "#E66101",
               "Site HS" = "#1B9E77",
               "Site AB" = "#B35806"), 
    name = "") +
  theme_void() +
  theme(panel.border = element_rect(color = "black", linewidth = 0.5, fill = NA),
    plot.margin = margin(0, 0, 0, 0, "mm"),
    legend.position = "none")

# Combine using grid
library(grid)
png("figures/site_maps.png", width = 6, height = 8, units = "in", res = 350)
print(main_map)
print(inset_map, vp = viewport(x = 0.67, y = 0.22, width = 0.25, height = 0.25))
dev.off()





