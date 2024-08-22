data {
  int<lower=0> N;
  vector [N] y;
  vector [N] x;
  
  real beta_mean_prior;
}

parameters {
  real intercept;
  real beta;
  real<lower=0> sigma;
}

model {
  y ~ normal(intercept + beta * x, sigma);
  intercept ~ normal(0,1);
  beta ~ normal(beta_mean_prior, 1);
  sigma ~ exponential(1);

}
