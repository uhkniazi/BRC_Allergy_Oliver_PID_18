**Diagnostic Model 2** utilised data collected for IgE sensitisation to the 6 peanut allergen components included on the ImmunoCAP ISAC and, as before, was built with data from both PA (n=50) and PS (n=40) individuals. A more robust version of the GLM previously described, was used in order to tolerate outliers . This is achieved by modifications to present the model as a finite mixture model, where one component of the mixture is a hierarchical logistic regression model, and the second component is a random ‘guessing’ parameter adjusting for outliers [Ref]. This provides a ‘weighted’ version of the estimated coefficients reducing the unwanted effects of outliers in the data.  

Kruschke JK. Doing Bayesian data analysis: A tutorial with R, JAGS, and Stan. 2nd ed; 2014
