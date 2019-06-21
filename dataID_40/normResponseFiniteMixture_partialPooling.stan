data {
  int<lower=1> Ntotal; // number of observations
  real y[Ntotal]; // response variable - normally distributed
  int<lower=1> Ncol; // total number of columns in model matrix
  matrix[Ntotal, Ncol] X; // model matrix
  int<lower=2> iMixtures; // number of mixture distributions
  //real iIntercepts[iMixtures];
}

parameters { // the parameters to track
  ordered[iMixtures] mu; // number of means to track Breaking the Labeling Degeneracy by Enforcing an Ordering (population Intercepts)
  real<lower=0.01> sigmaPop[iMixtures]; // scale parameters for normal distribution (population sigmas) 
  simplex[iMixtures] iMixWeights; // weights for the number of mixtures (should sum to one)
  // regression coefficients and other related parameters
  //real<lower=0.01> sigmaRan; // shared standard deviation for coefficients
  vector[Ncol] betas; // regression parameters
}
transformed parameters {
  vector[Ntotal] muFitted; // fitted value from linear predictor
  // fitted value
  muFitted = X * betas; 
  // assumption is that the effect sizes are the same in each distribution
  // see below for intercepts or population means
}
model {
  // see stan manual page 187 for an example
  real ps[iMixtures]; // temporary variable for log components
  //sigmaRan ~ cauchy(0, 2);
  betas ~ cauchy(0, 2);
  // any priors for mixture components go here 
  mu[1] ~ cauchy(0, 2);
  mu[2] ~ cauchy(0, 2);
  sigmaPop ~ cauchy(0, 2);
  iMixWeights ~ dirichlet(rep_vector(2.0, iMixtures));
  // loop to calculate likelihood
  for(n in 1:Ntotal){
    // second loop for number of mixture components
    for (k in 1:iMixtures){
      ps[k] = log(iMixWeights[k]) + normal_lpdf(y[n] | mu[k] + muFitted[n], sigmaPop[k]);
    }
    target += log_sum_exp(ps);
  }
}
