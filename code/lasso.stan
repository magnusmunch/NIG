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
  real<lower=0> lambdasq;
  vector<lower=0>[p] tausq;
  real<lower=0> sigmasq;
}

transformed parameters {
  vector<lower=0>[p] tau = sqrt(tausq);
  real<lower=0> lambda = sqrt(lambdasq);
  real<lower=0> sigma = sqrt(sigmasq);
}

model {
  // priors
  target += jeffreys_lpdf(sigmasq);   
  lambdasq ~ gamma(a0, b0);
  tausq ~ exponential(lambdasq/2);
  beta ~ normal(0., sigma*tau);
  
  // likelihood
  y ~ normal(x*beta, sigma);
}
