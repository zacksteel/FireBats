## Purpose: Organize data for BatAMP upload
## Project: FireBAts

library(tidyverse)
library(lubridate)

## Bring in count data
d <- read.csv("data/counts.csv") %>% 
  ## create some columns that match batamp template
  mutate(first_name = "Zack",
         last_name = "Steel",
         det_mfg = "SonoBat",
         det_model = "SM3BAT",
         ## Specify mic type
         mic_type = case_when(sm3_u1 == 0 ~ "SMM-U1",
                              sm3_u1 == 1 ~ "SM3-U1"),
         refl_type = "None",
         mic_ht_units = "meters",
         call_id_1 = "Sonobat 3",
         call_id_2 = NA,
         site_ID = forest,
         det_id = point,
         species = toupper(species),
         night = ymd(night),
         night = format(night, format="%m/%d/%Y")
  ) %>% 
  ## remove some columns
  select(-point, -sm3_u1, -forest)

## make it wide
d_w <- pivot_wider(d, names_from = species, values_from = ccounts) %>% 
  ##reorder columns
  select(first_name, last_name, y_coord, x_coord, det_mfg, det_model, mic_type, refl_type,
         mic_ht, mic_ht_units, everything())

## data must be within unique years
for(i in unique(d$year)) {
  sub <- filter(d_w, year == i) %>% 
    select(-year)
  write.csv(sub, paste0("data/batamp/zls_", i, ".csv"), row.names = F)
}
