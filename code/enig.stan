functions {
  // inverse Gaussian log probability density function
  real igauss_lpdf(real x, real theta, real lambda) {
    return -log(2*pi())/2 + log(lambda)/2 - 3*log(x)/2
           - lambda*square(x - theta)/(2*square(theta)*x);
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
  vector<lower=0>[p] ctalphainv;
  real<lower=0> lambda;
}

parameters {
  real beta0;
  vector[p] beta;
  real<lower=0> sigmasq;
  vector<lower=0>[p] gammasq;
}

transformed parameters {
  real<lower=0> sigma = sqrt(sigmasq);
  vector<lower=0>[p] gamma = sqrt(gammasq);
  vector<lower=0>[p] betasd = sigma*gamma;
  vector[n] mu;
  mu = beta0 + x*beta;
}

model {
  // priors
  target += jeffreys_lpdf(sigmasq);   
  for(j in 1:p) {
    beta[j] ~ normal(0., betasd[j]);
    target += igauss_lpdf(gammasq[j] | ctalphainv[j], lambda);
  }

  // likelihood
  y ~ normal(mu, sigma);
}
