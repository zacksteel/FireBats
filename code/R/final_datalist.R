## Purpose: Generate stan_data lists 
## Author: Zack Steel

library(tidyverse)

## Bring in looic table and filter for lowest value
lootab <- read.csv("data/lootab.csv") %>%
  group_by(sp) %>%
  filter(weight == max(weight)) %>%
  ungroup() %>%
  mutate(scale = as.factor(scale)) %>%
  as.data.frame()

## Read in data generating function
source("Code/R/model_data.R")

## list of stan_data objects
data_l <- lapply(1:nrow(lootab), function(x) 
  bat_stan_data(data_path = "data/occupancy_data.RData",
                sp = lootab[x,"sp"],
                sev_scale = lootab[x, "scale"]))
names(data_l) <- lootab$sp

## Save for later use
save(data_l, file = "data/stan_data.RData")

