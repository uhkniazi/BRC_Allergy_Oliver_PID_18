data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Ncol; // total number of columns in model matrix
  matrix[Ntotal, Ncol] X; // model matrix
  real y[Ntotal]; // response variable t distributed
  int<lower=0> Ncens; // number of censored data points 
  real<upper=min(y)>  rLower; // lower censoring point
  matrix[Ncens, Ncol] X2; // model matrix
}
// transformed data {
  // }
parameters {
  // parameters to estimate in the model
  vector[Ncol] betas; // regression parameters
  real<lower=0> sigmaPop; // population standard deviation
  real<lower=1> nu; // normality parameter for t distribution or degree of freedom 
  real<upper=rLower> y_cens[Ncens];
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  vector[Ncens] mu2; // fitted values from censored linear predictor
  // fitted value
  mu = X * betas; 
  mu2 = X2 * betas;
}
model {
  nu ~ exponential(1/29.0);
  // using diffuse prior
  sigmaPop ~ cauchy(0, 2);
  betas ~ cauchy(0, 2); //prior for the betas
  // likelihood function
  y ~ student_t(nu, mu, sigmaPop);
  y_cens ~ student_t(nu, mu2, sigmaPop);
}
generated quantities {
  vector[Ntotal] log_lik;
  for (i in 1:Ntotal) log_lik[i] = student_t_lpdf(y[i] | nu, mu[i], sigmaPop);
}