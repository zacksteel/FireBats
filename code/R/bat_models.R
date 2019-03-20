## Purpose: Run stan models for all species and scales
## Author: Zack Steel

library(tidyverse)
library(rstan)
library(rethinking)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) #allows for parallel computing

## Read in data generating function
source("Code/R/model_data.R")

## Set up species vector
spp <- c("anpa", "coto", "epfu", "labl", 
         "laci", "lano", "myca", "myci", "myev", "mylu",
         "myth", "myvo", "myyu", "pahe", 
         "eupe", "euma", "tabr") 

## Set up scale vector
sev_scales <- c(50, 100, 250, 500)


## Run occupancy models for each species and scale, saving after each species
#### Alternatively run just the selected models with stan_data from ..data/stan_data.RData
occ_mods <- list()

for(i in 1:length(spp)) {
  
  sp <- spp[i]
  
  cat("Running occupancy models for", sp, "\n")
  
  ## list of stan_data objects
  sd_l <- lapply(sev_scales, function(x) 
    bat_stan_data(data_path = "data/occupancy_data.RData",
                  sp = sp,
                  sev_scale = x))
  
  ## Compile model first to save time
  m <- stan_model("Code/stan/FireBatOcc.stan")
  
  ## For each severity scale, run a model
  mods_l <- lapply(1:length(sd_l), function(x) {
    
    stan_data <- sd_l[[x]]
    scale <- sev_scales[x]
    
    ## MCMC settings
    ni <- 2000
    nt <- 1
    nb <- 1000
    nc <- 3
    
    
    ## Parameters to (not) track
    pars <- c("b_site", "sigma_b_site", "sigma_b_fire", "tilde_b_site", "tilde_b_fire")
    
    ## Call Stan from R
    out <- sampling(m, data = stan_data, chains = nc, iter = ni,
                    pars = pars, include = F, seed = 1, 
                    control=list(adapt_delta=0.9),
                    open_progress = F)
    
    # =========================
    # = Timing and Efficiency =
    # =========================
    timing <- cbind(get_elapsed_time(out), rowSums(get_elapsed_time(out)))
    max_time <- max(timing)/60
    neff_md <- median(summary(out)$summary[,"n_eff"])
    rhat_md <- median(summary(out)$summary[,"Rhat"], na.rm=T) #returns NA for derived parameters like psi_con and z
    rhat_max <- max(summary(out)$summary[,"Rhat"], na.rm=T)
    eff <- neff_md/max_time
    
    time <- data.frame(tot_min = max_time, neff_md = neff_md, 
                       rhat_md = rhat_md, rhat_max = rhat_max, neff_min = eff)
    
    ## Return model and timing info
    list(model = out, scale = scale, data = stan_data, time = time)
  })
  
  names(mods_l) <- sev_scales
  
  ## Save and move on to next species
  occ_mods[[i]] <- mods_l
  names(occ_mods)[i] <- sp
  
  #### Modify file location as needed. Ultimately creates a file with models for 17 species x 4 scales. 
  #### File size is approximately 11 GB with current settings
  # save(occ_mods, file = "occMods.RData")
  
}

