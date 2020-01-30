functions {
  // Jeffrey's log probability density (up to constant)
  real jeffreys_lpdf(real x) {
    return -3*log(x)/2;
  }
}

data {
  int<lower=1> p;
  int<lower=1> n;
  matrix[n, p] x;
  vector[n] y; 
  real<lower=0> a0;
  real<lower=0> b0;
}

parameters {
  vector[p] beta;
  real<lower=0> lambda;
  real<lower=0> sigmasq;
}

transformed parameters {
  real<lower=0> sigma = sqrt(sigmasq);
  real<lower=0> betasd = sigma/sqrt(lambda);
}

model {
  // priors
  target += jeffreys_lpdf(sigmasq);   
  lambda ~ gamma(a0, b0);
  beta ~ normal(0., betasd);
  
  // likelihood
  y ~ normal(x*beta, sigma);
}
