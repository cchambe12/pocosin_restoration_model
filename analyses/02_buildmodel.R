#### Build Peat Model
## Started 25 October 2025 by Cat

### housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)

### Load Libraries
library(dplyr)
library(tidyr)
library(brms)
library(ggplot2)
library(ggsci)
library(bayesplot)
library(flextable)
library(webshot2)

### Set working directory
setwd("~/Documents/git/tnc_cprg/analyses/peatlands/")

set.seed(1221)

cleandat <- read.csv("output/clean_swails.csv") 


#### Evaluate collinearity
if(FALSE){
  plot(wtd ~ soiltemp, data = cleandat)
  plot(wtd ~ scPDSI, data = cleandat)
  plot(wtd ~ totalGDD, data = cleandat)
  plot(wtd ~ elevation, data = cleandat)
  plot(wtd ~ distkm, data = cleandat)
  plot(wtd ~ lai, data = cleandat)
  plot(wtd ~ ndvi, data = cleandat)
  
  plot(soiltemp ~ scPDSI, data = cleandat)
  plot(soiltemp ~ totalGDD, data = cleandat)
  plot(soiltemp ~ elevation, data = cleandat)
  plot(soiltemp ~ distkm, data = cleandat)
  plot(soiltemp ~ lai, data = cleandat)
  plot(soiltemp ~ ndvi, data = cleandat)
  
  plot(scPDSI ~ totalGDD, data = cleandat)
  plot(scPDSI ~ elevation, data = cleandat)
  plot(scPDSI ~ distkm, data = cleandat)
  plot(scPDSI ~ lai, data = cleandat)
  plot(scPDSI ~ ndvi, data = cleandat)
  
  plot(totalGDD ~ elevation, data = cleandat)
  plot(totalGDD ~ distkm, data = cleandat)
  plot(totalGDD ~ lai, data = cleandat)
  plot(totalGDD ~ ndvi, data = cleandat)
  
  plot(elevation ~ distkm, data = cleandat)
  plot(elevation ~ lai, data = cleandat)
  plot(elevation ~ ndvi, data = cleandat)
  
  plot(lai ~ distkm, data = cleandat)
  plot(lai ~ ndvi, data = cleandat)
  
  ggplot(aes(y = soilresp, x = wtd), data = cleandat) + geom_smooth() +
    geom_point() + theme_bw()
  
  ggplot(aes(y = soilresp, x = status), data = cleandat) + geom_smooth() +
    geom_point() + theme_bw() + geom_smooth(method="lm")
}


################################################################################
############################## Run the model ###################################
### Remove NAs
moddat <- cleandat %>%
  mutate(date = NULL,
         uniqueid = paste0(location, plot)) 
moddat <- moddat[complete.cases(moddat),]

### Adding in TX effect might capture elements missing from WTD and Soil Temp
### Simple model with random slopes for each site
soilrespmod = brm(soilresp ^ 0.25 ~ 
                    ### multiply variables that change after restoration
                    wtd.z + soiltemp.z + #   * tx.z
                    ### Add in individual site-level variables next
                    gdd.z +  lai.z + ndvi.z + 
                    pdsi.z + dist.z + elev.z +
                    ## Add in varying slopes
                    (wtd.z + soiltemp.z + 
                       gdd.z +  lai.z + ndvi.z + 
                       pdsi.z + dist.z + elev.z 
                     | uniqueid), 
                  data = moddat, control = list(adapt_delta=0.99, max_treedepth = 12),
                  chains = 4, cores = 2, iter = 3500, warmup = 2000
                  #prior = prior(exponential(1), class = "sd") ## weakly informative priors do not change results at all
                  ) 

soilrespmod.skew = brm(soilresp ~ 
                    ### multiply variables that change after restoration
                    wtd.z + soiltemp.z + #   * tx.z
                    ### Add in individual site-level variables next
                    gdd.z +  lai.z + ndvi.z + 
                    pdsi.z + dist.z + elev.z +
                    ## Add in varying slopes
                    (wtd.z + soiltemp.z + 
                       gdd.z +  lai.z + ndvi.z + 
                       pdsi.z + dist.z + elev.z 
                     | uniqueid), 
                  data = moddat, control = list(adapt_delta=0.99, max_treedepth = 12),
                  chains = 4, cores = 2, iter = 3500, warmup = 2000,
                  family = "skew_normal"
                  #prior = prior(exponential(1), class = "sd") ## weakly informative priors do not change results at all
) 

soilrespmod.noslopes = brm(soilresp ^ 0.25 ~ 
                    ### multiply variables that change after restoration
                    wtd.z + soiltemp.z + #   * tx.z
                    ### Add in individual site-level variables next
                    gdd.z +  lai.z + ndvi.z + 
                    pdsi.z + dist.z + elev.z +
                    ## Add in varying slopes
                    (1 | uniqueid), 
                  data = moddat, control = list(adapt_delta=0.99, max_treedepth = 12),
                  chains = 4, cores = 2, iter = 3500, warmup = 2000) 

### Simple model with random slopes for each site, without elevation
soilrespmod.noelev = brm(soilresp ^ 0.25 ~ 
                           ### multiply variables that change after restoration
                           wtd.z + soiltemp.z + #   * tx.z
                           ### Add in individual site-level variables next
                           gdd.z +  lai.z + ndvi.z + 
                           pdsi.z + dist.z + #elev.z +
                           ## Add in varying slopes
                           (wtd.z + soiltemp.z + 
                              gdd.z +  lai.z + ndvi.z + 
                              pdsi.z + dist.z #+ elev.z 
                            | uniqueid), 
                         data = moddat, control = list(adapt_delta=0.99, max_treedepth = 12),
                         chains = 4, cores = 2, iter = 3500, warmup = 2000) 

soilresp.simple = brm(soilresp ^ 0.25 ~ wtd.z + soiltemp.z,
                      data = moddat, control = list(adapt_delta=0.99, max_treedepth = 12),
                      chains = 4, cores = 2, iter = 3500, warmup = 2000)

soilresp.simple.wtd = brm(soilresp ^ 0.25 ~ wtd.z,
                      data = moddat, control = list(adapt_delta=0.99, max_treedepth = 12),
                      chains = 4, cores = 2, iter = 3500, warmup = 2000)

### Save the model
save(soilrespmod, file = "models/soilresp.Rdata")
save(soilrespmod.noslopes, file = "models/soilresp.noslopes.Rdata")
save(soilrespmod.noelev, file = "models/soilresp.noelev.Rdata")
save(soilresp.simple, file = "models/soilresp.simple.Rdata")
save(soilresp.simple.wtd, file = "models/soilresp.simple.wtd.Rdata")


if(FALSE){ ## Complete model checks
  
  load("models/soilresp.Rdata")
  load("models/soilresp.noslopes.Rdata")
  load("models/soilresp.noelev.Rdata")
  load("models/soilresp.simple.Rdata")
  load("models/soilresp.simple.wtd.Rdata")
  
  ### Evaluate model performance
  # Extract neff_ratio for fixed effects only
  compare_neff <- function(model, model_name) {
    pars <- c("sigma", grep("^b_", variables(model), value = TRUE))
    neff <- neff_ratio(model, pars = pars)
    
    data.frame(
      variable = names(neff),
      neff_ratio = as.numeric(neff),
      model = model_name
    )
  }
  
  # Get data for each model
  df_full <- compare_neff(soilrespmod, "full")
  df_noslopes <- compare_neff(soilrespmod, "noslopes")
  df_noelev <- compare_neff(soilrespmod.noelev, "noelev")
  df_simple <- compare_neff(soilresp.simple, "simple")
  df_wtd <- compare_neff(soilresp.simple.wtd, "wtd")
  
  # Pivot to side-by-side format
  comparison <- merge(df_full, df_noslopes, by = "variable", all = TRUE)
  comparison <- merge(comparison, df_noelev, by = "variable", all = TRUE)
  comparison <- merge(comparison, df_simple, by = "variable", all = TRUE)
  comparison <- merge(comparison, df_wtd, by = "variable", all = TRUE)
  
  # Clean up column names
  names(comparison) <- c("variable", "full", "model.x", "noslopes", "model.x1",
                         "noelev", "model.y", 
                         "simple", "model.z", "wtd", "model")
  comparison <- comparison[, c("variable", "full", "noslopes", "noelev", "simple", "wtd")]
  comparison[, 2:6] <- round(comparison[, 2:6], 2)
  
  comparison %>%
    flextable() %>%
    set_header_labels(variable = "Parameter", full = "Updated", noslopes = "No Slopes",
                      noelev = "No Elevation", gamma = "Gamma Distribution",
                      simple = "Baseline", wtd = "WTD Only") %>%
    theme_booktabs() %>%
    autofit() %>%
    save_as_image(path = "figures/neff_modelcompare.png", width = 7, height = 6)
  
  # Comprehensive summary
  summary(soilrespmod)
  summary(soilrespmod.skew)
  summary(soilrespmod.noslopes) ## slighly higher sigma
  summary(soilrespmod.noelev) ## has one parameter with Rhat of 1.02
  summary(soilresp.simple)
  summary(soilresp.simple.wtd)
  
  # Compare posterior predictive performance
  p_all <- ppc_dens_overlay(y = moddat$soilresp, yrep = posterior_predict(soilrespmod, draws = 50)) +
    ggtitle("Full Model")
  
  p_skew <- ppc_dens_overlay(y = moddat$soilresp, yrep = posterior_predict(soilrespmod.skew, draws = 50)) +
                               ggtitle("Skew Model")
  
  p_noslopes <- ppc_dens_overlay(y = moddat$soilresp, yrep = posterior_predict(soilrespmod.noslopes, draws = 50)) +
    ggtitle("Full Model - no slopes")
  
  p_noelev <- ppc_dens_overlay(y = moddat$soilresp, yrep = posterior_predict(soilrespmod.noelev, draws = 50)) +
    ggtitle("Full Model - no elev")
  
  p_simp <- ppc_dens_overlay(y = moddat$soilresp, yrep = posterior_predict(soilresp.simple, draws = 50)) +
    ggtitle("Simple Model (WTD + SoilTemp)")
  
  p_simp.wtd <- ppc_dens_overlay(y = moddat$soilresp, yrep = posterior_predict(soilresp.simple.wtd, draws = 50)) +
    ggtitle("WTD Only")
  
  gridExtra::grid.arrange(p_all, p_noslopes, p_noelev, p_simp, p_simp.wtd, nrow = 1)
  
  # Add to your libraries
  library(loo)
  
  # Calculate LOO for each model
  loo1 <- add_criterion(soilrespmod, "loo")
  loo1a <- add_criterion(soilrespmod.noslopes, "loo")
  loo2 <- add_criterion(soilrespmod.noelev, "loo")
  loo3 <- add_criterion(soilresp.simple, "loo")
  loo4 <- add_criterion(soilresp.simple.wtd, "loo")
  
  # Compare them
  loo_comparison <- loo_compare(loo1, loo1a, loo2, loo3, loo4, criterion = "loo")
  
  # Convert to data frame and rename
  loo_df <- as.data.frame(loo_comparison)
  loo_df <- loo_df %>%
    mutate(loos = rownames(loo_df),
           model = case_when(
             loos == "loo1" ~ "Updated model",
             loos == "loo1a" ~ "Updated model - no slopes",
             loos == "loo2" ~ "Updated model - no elevation",
             loos == "loo3" ~ "Baseline model",
             loos == "loo4" ~ "Baseline model - WTD only"))
       
  loo_df <- loo_df[, c("model", "elpd_diff", "se_diff")]
  
  # Round to 2 decimal places
  loo_df$elpd_diff <- round(loo_df$elpd_diff, 2)
  loo_df$se_diff <- round(loo_df$se_diff, 2)
  
  # Rename columns for clarity
  colnames(loo_df) <- c("Model", "elpd_diff", "SE_diff")
  
  loo_df %>%
    flextable() %>%
    theme_booktabs() %>%
    autofit() %>%
    save_as_image(path = "figures/loo_modelcompare.png", width = 7, height = 6)
  
  ### Check out pareto-k more specifically
  plot(loo(soilrespmod), diagnostic = c("k", "ESS", "n_eff"),
    label_points = FALSE, main = "PSIS diagnostic - Main Model")
  plot(loo(soilrespmod.skew), diagnostic = c("k", "ESS", "n_eff"),
       label_points = FALSE, main = "PSIS diagnostic - Skewed")
  
  ## Show trace plots for each
  plot(soilrespmod) 
  plot(soilrespmod.noslopes)
  plot(soilrespmod.noelev)
  plot(soilresp.simple)
  plot(soilrespmod.simple.wtd) 

# In-sample: R-squared
r2_full <- bayes_R2(soilrespmod)
r2_full.skew <- bayes_R2(soilrespmod.skew)
r2_full_noslopes <- bayes_R2(soilrespmod.noslopes)
r2_noelev <- bayes_R2(soilrespmod.noelev)
r2_simple <- bayes_R2(soilresp.simple)
r2_wtd <- bayes_R2(soilresp.simple.wtd)

r2_full; r2_full_noslopes; r2_noelev; r2_simple; r2_wtd

r2_full %>% as.data.frame() %>% 
  mutate(Model = "Updated model") %>% select(Model, Estimate) %>%
  bind_rows(r2_full_noslopes %>% as.data.frame() %>%
              mutate(Model = "Updated model - no slopes") %>% select(Model, Estimate),
            r2_noelev %>% as.data.frame() %>%
              mutate(Model = "Updated model - no elevation") %>% select(Model, Estimate),
            r2_simple %>% as.data.frame() %>%
              mutate(Model = "Baseline model") %>% select(Model, Estimate),
            r2_wtd %>% as.data.frame() %>%
              mutate(Model = "Baseline model - WTD only") %>% select(Model, Estimate)) %>%
  mutate(R2 = round(Estimate, digits = 2)) %>%
  select(Model, R2) %>%
  flextable() %>%
  theme_booktabs() %>%
  autofit() %>%
  save_as_image(path = "figures/R2_modelcompare.png", width = 7, height = 6)



pp_check(soilrespmod, stat = "mean", type = "stat") ## min, max, and sd look good too
pp_check(soilrespmod.skew, stat = "sd", type = "stat") ## all look bad
pp_check(soilrespmod.noslopes, stat = "mean", type = "stat")
pp_check(soilrespmod.noelev, stat = "mean", type = "stat")
pp_check(soilresp.simple, stat = "mean", type = "stat")
pp_check(soilresp.simple.wtd, stat = "mean", type = "stat")

}















