data {
  int<lower=0> N;
  vector [N] detect_index;
  vector [N] varprop_index;
}

parameters {
  real intercept;
  real BETA;
  real<lower=0> sigma;
}

model {

    varprop_index ~ normal(intercept + BETA * detect_index, sigma);
  
  
  intercept ~ normal(0,1);
  BETA ~ normal(1, 1);
  sigma ~ exponential(1);

}
