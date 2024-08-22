data {
  int<lower=0> N;
  vector [N] route_index;
  vector [N] point_index;
}

parameters {
  real intercept;
  real beta;
  real<lower=0> sigma;
}

model {
  route_index ~ normal(intercept + beta * point_index, sigma);
  intercept ~ normal(0,1);
  beta ~ normal(50, 1);
  sigma ~ exponential(1);

}
