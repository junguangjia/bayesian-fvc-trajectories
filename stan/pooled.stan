data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  vector[N] y;
  vector[K] beta_prior_mean;
  vector<lower=0>[K] beta_prior_sd;
  real<lower=0> sigma_prior_sd;
}
parameters {
  vector[K] beta;
  real<lower=0> sigma;
}
model {
  beta ~ normal(beta_prior_mean, beta_prior_sd);
  sigma ~ normal(0, sigma_prior_sd);
  y ~ normal(X * beta, sigma);
}
