// This code was developed by David Bauman (david.bauman@oxfordecosystems.co.uk)
// last version in 25.07.2024

// Based on https://vasishth.github.io/bayescogsci/book/ch-complexstan.html 
// Here we use a matrix of parameters for intercept and predictors. This will be the basis to then include the covariance matrix.

// For this vectorised specification of the model to work:
// - The predictors need to be vector[n] pred (and not real pred[n]);
//

// Notes regarding addition of a covariance matrix:
// - When declaring z_u_sp, we use [3, n_sp] instead of [n_sp, 3] used before (because will be transposed, in transformed parameters block); 
// - but for u_sp, we still use [n_sp, 3] (because z_u_sp was transposed)


data {
int<lower=0> n;         // nb observations
int<lower=0> n_plot;    // nb of plots
int<lower=0> n_sp;      // nb of species
//int<lower=0> n_interval;      // nb of time intervals b/ censuses
int mort[n];            // death status (1, dead, 0, survived)
vector[n] dbh;            // dbhs of observations
// real dbh[n];           // dbhs of observations
vector[n] gr;             // grs of observations
// real gr[n];             // growth rates
vector[n] time;           // time interval, in years
// real time[n];           // time interval, in years
// array[n] int<lower = 1, upper = n_sp> sp; // species ID of the n observations
int sp[n];           // species ID of the n observations
// array[n] int<lower = 1, upper = n_plot> plot; // plot ID of the n observations
int plot[n];         // plot id for the corresponding n obs.
// array[n] int<lower = 1, upper = n_interval> interval; // interval ID of the n observations
//int interval[n];         // interval id for the corresponding n obs.
}

parameters {
  // Group-level parameters:
  matrix[3, n_sp] z_u_sp;      // matrix for the z-score sp-level params; 3 columns for intercept (1) + nb of predictors (2, here: "dbh", then "gr") 
  // vector[n_sp] a_sp_z;      // "_z" as in "z-score transformed", to be used to recover the parameter on scale of interest, in transf. params.
  // vector[n_sp] b_dbh_sp_z;
  // vector[n_sp] b_gr_sp_z;
  cholesky_factor_corr[3] L_u_sp; 

  vector[n_plot] a_plot_z; 
  //vector[n_interval] a_interval_z; 
  
  // Hyperparameters:
  real a;
  real b_dbh;
  real b_gr;
  vector<lower = 0>[3] sigma_sp;  // vector for the sp-level SD params (3 for intercept + the 2 predictors)
  // real<lower=0> sigma_a_sp;      // for sp-level varying intercept
  // real<lower=0> sigma_b_dbh;     // for sp-level varying slopes
  // real<lower=0> sigma_b_gr;      // for sp-level varying slopes

  real<lower=0> sigma_a_plot;    // for plot-level varying intercept
  // real<lower=0> sigma_a_interval; // for interval-level varying intercept
}

transformed parameters {
  // Group-level parameters:
  matrix[n_sp, 3] u_sp;  // matrix for the back-transformed sp-level params
  vector[n_plot] a_plot; 
  //vector[n_interval] a_interval; 
  
  // Non-centered priors (smuggle the hyperparameters out of the adaptive priors):
  u_sp = (diag_pre_multiply(sigma_sp, L_u_sp) * z_u_sp)'; // note the transpose here
  a_plot = a_plot_z * sigma_a_plot;
  // a_interval = a_interval_z * sigma_a_interval;
}

model {
  
  // Declare probability of death 'p':
  real p_logit; // p on logit scale
  real p;       // p on natural scale (0 to 1)
  
  // Non-centered priors, to avoid having hyperparameters in the adaptive priors (often more efficient MCMC sampling when limited sample size):
  target += std_normal_lpdf(to_vector(z_u_sp));
  target += lkj_corr_cholesky_lpdf(L_u_sp | 2);

  target += std_normal_lpdf(a_plot_z);
  // target += std_normal_lpdf(a_interval_z);
  // a_plot_z ~ normal(0, 1);
  // a_interval_z ~ normal(0, 1);
  
  // Priors on hyperparameters for varying effects:
  target += normal_lpdf(a | -1, 2);  // hyperprior on the grand mean in a_sp ~ normal(a, sigma_a_sp), on logit scale: -3 corresponds to a mort. prob. of 0.04 and 1 to 0.73

  target += normal_lpdf(b_dbh | 0, 1);  
  target += normal_lpdf(b_gr | 0, 1);  
  // b_dbh ~ normal(0, 1);
  // b_gr ~ normal(0, 1);
  
  target += normal_lpdf(sigma_sp | 0, 1) - 3 * normal_lccdf(0 | 0, 1); // we multiply by 3 because there are three PDFs that need to be corrected for the truncation
  // Equivalent to:
  // target += normal_lpdf(sigma_sp[1] | 0, 1) - normal_lccdf(0 | 0, 1); 
  // target += normal_lpdf(sigma_sp[2] | 0, 1) - normal_lccdf(0 | 0, 1); 
  // target += normal_lpdf(sigma_sp[3] | 0, 1) - normal_lccdf(0 | 0, 1); 

  target += normal_lpdf(sigma_a_plot | 0, 1) - normal_lccdf(0 | 0, 1);
  //target += normal_lpdf(sigma_a_interval | 0, 1) - normal_lccdf(0 | 0, 1);

  // Likelihood:
    for (i in 1:n) {
    p_logit = a + u_sp[sp[i], 1] + (b_dbh + u_sp[sp[i], 2]) .* dbh[i] + (b_gr + u_sp[sp[i], 3]) .* gr[i] + 
                a_plot[plot[i]];// + a_interval[interval[i]];
    p = inv_logit(p_logit);
    // surv[i] ~ bernoulli(p^time[i]); // When modelling survival (surv: 1 = survided; 0 = died); Where p is the survival probability
    target += bernoulli_lpmf(mort[i] | 1-((1-p)^time[i]));  // When modelling mortality (mort: 1 = died; 0 = survived); Where p is the death probability
  }
}

generated quantities {
  real rho_a_dbh = (L_u_sp * L_u_sp')[1, 2];
  real rho_a_gr = (L_u_sp * L_u_sp')[1, 3];
  real rho_dbh_gr = (L_u_sp * L_u_sp')[2, 3];
  // The above three correlations could be extracted in one line, defining a matrix rho instead of three real:
  // matrix[2, 2] Omega;
  // Omega = multiply_lower_tri_self_transpose(L_u_sp);
  // Generate posteriors of species-level responses to dbh and gr (on the link-function scale), summing hyperparam. and sp deviations:
  vector[n_sp] effect_dbh_by_sp = b_dbh + u_sp[ , 2];
  vector[n_sp] effect_gr_by_sp = b_gr + u_sp[ , 3];
}
