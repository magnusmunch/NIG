functions {
  // inverse Gaussian log probability density function
  real igauss_lpdf(real x, real theta, real lambda) {
    return -log(2*pi())/2 + log(lambda)/2 - 3*log(x)/2
           - lambda*square(x - theta)/(2*theta^2*x);
  }
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
  real<lower=0> theta;
  real<lower=0> lambda;
}

parameters {
  real beta0;
  vector[p] beta;
  real<lower=0> sigmasq;
  real<lower=0> gammasq;
}

transformed parameters {
  real sigma = sqrt(sigmasq);
  real gamma = sqrt(gammasq);
  real betavar = sigma*gamma;
  vector[n] mu;
  mu = beta0 + x*beta;
}

model {
  // priors
  target += jeffreys_lpdf(sigmasq);   
  target += igauss_lpdf(gammasq | theta, lambda);
  beta ~ normal(0., betavar);
  // likelihood
  y ~ normal(mu, sigma);
}
