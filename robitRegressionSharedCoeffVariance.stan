data {
    int<lower=1> Ntotal; // number of observations
    int<lower=1> Ncol; // total number of columns in model matrix
    matrix[Ntotal, Ncol] X; // model matrix
    int y[Ntotal]; // response variable binomial distributed
}
// transformed data {
// }
parameters {
  // parameters to estimate in the model
    vector[Ncol] betas; // regression parameters
    //vector<lower=0>[(Ncol-1)] sigmas; // scale parameters
    real<lower=0> tau; // standard deviation for deflections
    //real<lower=3> nu;  // degrees of freedom for t-distribution
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  vector[Ncol] betas2;
  betas2[1] = betas[1];
  // tau can be close to zero henece parameterized this way
  betas2[2:Ncol] = betas[2:Ncol] * tau; 
  mu = X * betas2; 
}
model {
  vector[Ntotal] pi; 
  //sigmas ~ cauchy(0, 1);
  tau ~ exponential(1);
  betas[1] ~ cauchy(0, 10); //prior for the betas
  betas[2:Ncol] ~ normal(0, 1); // 
  //nu ~ uniform(3, 100);
  
  // latent variable formulation
  for (i in 1:Ntotal){
    pi[i] = student_t_cdf(mu[i], 4, 0, 1);
    // likelihood function
    y[i] ~ bernoulli(pi[i]);
  } 
}
