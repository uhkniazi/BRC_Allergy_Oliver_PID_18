data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Ncol; // total number of columns in model matrix
  matrix[Ntotal, Ncol] X; // model matrix
  real y[Ntotal]; // response variable t distributed
}
// transformed data {
  // }
parameters {
  // parameters to estimate in the model
  vector[Ncol] betas; // regression parameters
  real<lower=0> sigmaPop; // population standard deviation
  real<lower=0> sigmaRan; // shared standard deviation for coefficients
  real<lower=1> nu; // normality parameter for t distribution or degree of freedom 
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  // fitted value
  mu = X * betas; 
}
model {
  nu ~ exponential(1/29.0);
  // using diffuse prior
  sigmaPop ~ cauchy(0, 2);
  sigmaRan ~ cauchy(0, 2);
  betas[1] ~ cauchy(0, 2); //prior for the intercept
  betas[2:Ncol] ~ normal(0, sigmaRan);
  // likelihood function
  y ~ student_t(nu, mu, sigmaPop);
}
