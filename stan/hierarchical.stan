data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> K;
  array[N] int<lower=1, upper=J> patient;
  matrix[N, K] X;
  vector[N] time;
  vector[N] y;
  vector[K] beta_prior_mean;
  vector<lower=0>[K] beta_prior_sd;
  vector<lower=0>[2] tau_prior_sd;
  real<lower=0> sigma_prior_sd;
}
parameters {
  vector[K] beta;
  vector<lower=0>[2] tau;
  cholesky_factor_corr[2] L_Omega;
  matrix[2, J] z;
  real<lower=0> sigma;
}
transformed parameters {
  matrix[2, J] u = diag_pre_multiply(tau, L_Omega) * z;
}
model {
  vector[N] mu = X * beta;
  beta ~ normal(beta_prior_mean, beta_prior_sd);
  tau ~ normal(0, tau_prior_sd);
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(z) ~ std_normal();
  sigma ~ normal(0, sigma_prior_sd);
  for (n in 1:N) {
    mu[n] += u[1, patient[n]] + u[2, patient[n]] * time[n];
  }
  y ~ normal(mu, sigma);
}
generated quantities {
  // In a 2 x 2 Cholesky correlation factor, its [2, 1] entry is rho.
  real rho = L_Omega[2, 1];
}
