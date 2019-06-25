data {
  int<lower=1> Ntotal; // number of observations
  int<lower=1> Ncol; // total number of columns in model matrix
  int<lower=1> NscaleBatches;
  matrix[Ntotal, Ncol] X; // model matrix
  real y[Ntotal]; // response variable t distributed
  int<lower=1, upper=NscaleBatches> NBatchMap[(Ncol-1)]; // mapping variable to model matrix columns
  int<lower=0> Ncens; // number of censored data points 
  real<upper=min(y)>  rLower; // lower censoring point
  matrix[Ncens, Ncol] X2; // model matrix
}
// transformed data {
  // }
parameters {
  // parameters to estimate in the model
  vector[Ncol] betas; // regression parameters
  real<lower=0> sigmaPop; // population level scale for t distribution
  real<lower=0.01> sigmaRan[NscaleBatches]; // shared standard deviation for coefficients
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
  real sigmaRan_expanded[(Ncol-1)];
  nu ~ exponential(1/29.0);
  // using diffuse prior
  sigmaPop ~ cauchy(0, 2);
  sigmaRan ~ cauchy(0, 2);
  betas[1] ~ cauchy(0, 2); //prior for the intercept
  // vector expansion by mapping to a larger vector/array
  sigmaRan_expanded = sigmaRan[NBatchMap];
  betas[2:Ncol] ~ normal(0, sigmaRan_expanded);
  // likelihood function
  y ~ student_t(nu, mu, sigmaPop);
  y_cens ~ student_t(nu, mu2, sigmaPop);
}
