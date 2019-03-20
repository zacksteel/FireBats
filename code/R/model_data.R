## Purpose: Set up Stan_data object for postfire bat occupancy 
## Author: Zack Steel

bat_stan_data <- function(data_path, # path to preppred RData file
                          sp,        # species code
                          sev_scale) # severity scale
  {
  
  library(tidyverse)
  library(rstan)
  library(rethinking)
  
  ## Bring in prepped data from OccMOdData.R
  load(data_path)
  
  ## Change observation data to something generic
  obs_d <- obs_bin
  obs_d <- mutate(obs_d, resp = present) %>%
    select(-present)
  rm(obs_bin)
  
  ## scaling function
  f_s <- function(x) {(x - mean(x, na.rm=T)) / sd(x, na.rm=T)}
  
  ## subset and mutate
  d <- filter(obs_d, species == sp) %>%  
    arrange(point, year, visit)
  sc <- filter(sev, scale == sev_scale) %>%
    merge(site_covs, by = c("point", "fire")) %>%
    mutate(point = factor(point)) %>%
    ## Standardize all covariates & and add quadratics where appropriate
    mutate(sev_s = f_s(sevmn),
           sevsq_s = sev_s * sev_s,
           sevsdreal_s = f_s(sevsd),
           sevsd_s = f_s(sd_resid), 
           distw_s = f_s(distw),
           dwlog_s = f_s(dw_log),
           elev_s = f_s(elev),
           cc_s = f_s(cc15),
           newmic = ifelse(sm3_u1 == 1, 0, 1),
           mic = as.integer(newmic), 
           ## Also create an integer category for fire, site 
           fireid = fire,
           fire = as.integer(fire),
           site = as.integer(point)) %>%
    select(point, year, site, fireid, fire, 
           sevmn, sevsd, distw, dw_log, elev, cc15, 
           sev_s, sevsq_s, sevsd_s,  
           distw_s, dwlog_s, elev_s, cc_s, mic) %>% 
    arrange(point, year) %>%
    as.data.frame()
  
  oc <- mutate(obs_covs,
               jday_s = f_s(jday),
               jday_s2 = jday_s^2,
               st_s = f_s(scr_tri),
               temp_s = f_s(temp_avg)) %>%
    arrange(point, year)
  
  ## make obs and obs covs wide with M sites as rows and J visits as columns
  ## Spread observation covariates for multiple surveys
  ## have to do this for one variable at a time
  jday <- select(oc, fire, point, year, visit, jday_s) %>% 
    dplyr::rename(jday = visit) %>%
    ## spread by visit number
    spread(key = jday, value = jday_s, sep = "_") %>%
    arrange(point, year)
  
  jday2 <- select(oc, fire, point, year, visit, jday_s2) %>%
    dplyr::rename(jday2 = visit) %>%
    ## spread by visit number
    spread(key = jday2, value = jday_s2, sep = "_") %>%
    arrange(point, year)
  
  noise <- select(oc, fire, point, year, visit, st_s) %>%
    dplyr::rename(noise = visit) %>%
    ## spread by visit number
    spread(key = noise, value = st_s, sep = "_") %>%
    arrange(point, year)
  
  temp_mn <- select(oc, fire, point, year, visit, temp_s) %>%
    dplyr::rename(temp_mn = visit) %>%
    ## spread by visit number
    spread(key = temp_mn, value = temp_s, sep = "_") %>%
    arrange(point, year)
  
  ## Spread species observations as well
  ## Activity
  d_w <- select(d, -night) %>%
    spread(key = visit, value = resp, sep = "_") 
  
  ## make sure everything still lines up
  all.equal(d_w[,c("point","year")], sc[,c("point","year")])
  all.equal(d_w[,c("point","year")], jday[,c("point","year")])
  all.equal(d_w[,c("point","year")], noise[,c("point","year")])
  all.equal(d_w[,c("point","year")], temp_mn[,c("point","year")])
  
  ## Just the obs data 
  y <- select(d_w, -point, -year, -fire, -species) %>%
    as.matrix()
  jday_oc <- select(jday, -point, -year, -fire) %>%
    as.matrix()
  jday2_oc <- select(jday2, -point, -year, -fire) %>%
    as.matrix()
  noise_oc <- select(noise, -point, -year, -fire) %>%
    as.matrix()
  temp_oc <- select(temp_mn, -point, -year, -fire) %>%
    as.matrix()
  
  ## Create an indicator vector of when visits stopped
  last <- sapply(1:nrow(y), #same for all species
                 function(i) max(grep(FALSE, is.na(y[i,]))))
  ## Replace NAs with zeros that will be ignored according to 'last'
  y[is.na(y)] <- 0
  jday_oc[is.na(jday_oc)] <- 0
  jday2_oc[is.na(jday2_oc)] <- 0
  noise_oc[is.na(noise_oc)] <- 0
  temp_oc[is.na(temp_oc)] <- 0
  
  ## Bundle and summarize stan data
  stan_data <- list(y = y, #Observations
                    ## Random variables
                    site = sc$site,
                    fire = sc$fire,
                    ## Continuous covariates on occupancy
                    sev = sc$sev_s,
                    sevsq = sc$sevsq_s,
                    sevsd = sc$sevsd_s,
                    elev = sc$elev_s,
                    distw = sc$dwlog_s,
                    ## Covariates on detection
                    cc = sc$cc_s,
                    mic = sc$mic,
                    jday = jday_oc,
                    jday2 = jday2_oc,
                    noise = noise_oc,
                    temp = temp_oc,
                    ## Lengths for indexing
                    M = dim(y)[1],
                    J = dim(y)[2],
                    N_site = max(sc$site),
                    N_fire = max(sc$fire),
                    last = last,
                    ## Also add original data
                    sc = sc,
                    oc = oc,
                    obs_w = d_w)
  
}