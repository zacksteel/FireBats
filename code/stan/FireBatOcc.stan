// Purpose: Site-occupancy model for effects of burn severity on bats
// Author: Zack Steel

data {
  int<lower=1> M;                 // Number of sites x years
  int<lower=1> N_site;            // Number of unique sites
  int<lower=1> N_fire;            // Number of fires
  int<lower=1> J;                 // Number of temporal replications
  int<lower=0,upper=1> y[M, J];   // Observations
  
  int<lower=1,upper=N_site> site[M];                 // Covariates
  int<lower=1,upper=N_fire> fire[M];
  vector[M] sev;                  // Mean basal area mortality
  vector[M] sevsq;                // sev squared
  vector[M] sevsd;                // pyrodiversity
  vector[M] elev;                 // elevation
  vector[M] distw;                // distance to water
  vector[M] cc;                   // canopy cover
  vector[M] mic;                  // indicator of SMM_U1 microphone
  
  matrix[M, J] jday;              // Julian day
  matrix[M, J] jday2;             // Julian day squared
  matrix[M, J] noise;             // proportion of files recorded that were noise
  matrix[M, J] temp;              // mean temperature
  
  int last[M];                    // indicator of how many visits for each site
}

transformed data {
  int<lower=0,upper=J> sum_y[M];  // Number of occupation for each site
  int<lower=0,upper=M> occ_obs;   // Number of observed occupied sites

  occ_obs = 0;
  for (i in 1:M) {
    sum_y[i] = sum(y[i]);
    if (sum_y[i])
      occ_obs = occ_obs + 1;
  }
}

parameters {
  real beta;
  real<lower=0> sigma_b_site;
  real<lower=0> sigma_b_fire;
  vector[N_site] tilde_b_site;
  vector[N_fire] tilde_b_fire;
  
  real b_distw;
  real b_elev;
  real b_sev;
  real b_sevsq;
  real b_sevsd;
  
  real alpha;
  
  real a_jday;
  real a_jday2;
  real a_cc;
  real a_mic;
  real a_noise;
  real a_temp;
}

transformed parameters {
  // non-centered intercepts
  vector[N_site] b_site = sigma_b_site * tilde_b_site;
  vector[N_fire] b_fire = sigma_b_fire * tilde_b_fire;
  
  vector[M] logit_psi;            // Logit occupancy probability
  matrix[M, J] logit_p;           // Logit detection probability
  vector[M] log_lik;              // vector to store log likelihood
  
  for (i in 1:M) {
    logit_psi[i] = beta + b_site[site[i]] + b_fire[fire[i]] +// global intercept & random intercepts
      b_distw * distw[i] + b_elev * elev[i] + // parameters estimated with all data
      b_sev * sev[i] + b_sevsd * sevsd[i] + b_sevsq * sevsq[i];
  }

      
  for (i in 1:M) {
    for (j in 1:J) {
      logit_p[i,j] = alpha + 
        a_jday * jday[i,j] + a_jday2 * jday2[i, j] + 
        a_cc * cc[i] + a_mic * mic[i] +
        a_noise * noise[i,j] + a_temp * temp[i,j];
    }
  }
  
  for (i in 1:M) {
    if (sum_y[i]) {
      // observed, present
      log_lik[i] = bernoulli_logit_lpmf(1 |  logit_psi[i])
                   + bernoulli_logit_lpmf(y[i, 1:last[i]] | logit_p[i, 1:last[i]]);
    } else {
      // not observed
      log_lik[i] = log_sum_exp( // either:
                     // present and not detected
                     bernoulli_logit_lpmf(1 | logit_psi[i])
                       + bernoulli_logit_lpmf(0 | logit_p[i, 1:last[i]]),
                     // or absent
                     bernoulli_logit_lpmf(0 | logit_psi[i]));
    }
  }
}



model {
  
  // Priors
  beta ~ normal(0, 1.5);
  
  sigma_b_fire ~ normal(0, 1);
  tilde_b_fire ~ normal(0, 1);
  
  sigma_b_site ~ normal(0, 1);
  tilde_b_site ~ normal(0, 1);
  
  b_distw ~ normal(0, 1);
  b_elev ~ normal(0, 1);
  b_sev ~ normal(0, 1);
  b_sevsq ~ normal(0, 1);
  b_sevsd ~ normal(0, 1);
  
  alpha ~ normal(0, 1);
  
  a_jday ~ normal(0, 1);
  a_jday2 ~ normal(0, 1);
  a_cc ~ normal(0, 1);
  a_mic ~ normal(0, 1);
  a_noise ~ normal(0, 1);
  a_temp ~ normal(0, 1);
  
  // increment the log density
  target += sum(log_lik);

}

generated quantities {
  int occ_fs;       // Number of occupied sites
  real psi_con[M];  // prob occupied conditional on data
  int z[M];         // occupancy indicator, 0/1
  
  for (i in 1:M) {
    if (sum_y[i] == 0) {  // species not detected
      real psi = inv_logit(logit_psi[i]);
      vector[J] q = inv_logit(-logit_p[i])';  // q = 1 - p
      real qT = prod(q[]);
      psi_con[i] = (psi * qT) / (psi * qT + (1 - psi));
      z[i] = bernoulli_rng(psi_con[i]);
    } else {             // species detected at least once
      psi_con[i] = 1;
      z[i] = 1;
    }
  }
  occ_fs = sum(z);
}
