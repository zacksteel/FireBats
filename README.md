# FireBats
Code accompanying 2019 manuscript

R functions in code/R include:  

- *bat_models.R* which can be used to run occupancy models for all 17 species and 4 spatial scales  

- *final_datalist* compiles stan_data lists for each species and the scale selected by LOOIC  
- *model_data.R* the two functions above depend on model_data.R to format bat occurrence and environmental data for stan models  

*code/stan/FireBatOcc.stan* contains the occupancy model coded in Stan. This is run in bat_models.R  

Data files include:  

- *occupancy_data.RData* contains all the bat occurrence and environmental data used in the manuscript analysis  

- *stan_data.RData* contains the data for each species and scale selected by LOOIC, formated for the Stan model  

- *lootab.csv* is the looic table used for model selection  

- *OccPars.csv* is the parameter estimates from the final model used in the manuscript  
