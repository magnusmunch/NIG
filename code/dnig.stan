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
  real<lower=0> ctalphainv;
  real<lower=0> lambda;
  real<lower=0> ztthetainv;
  real<lower=0> kappa;
}

parameters {
  real beta0;
  vector[p] beta;
  real<lower=0> sigmasq;
  real<lower=0> gammasq;
  vector<lower=0>[p] tausq;
}

transformed parameters {
  real<lower=0> sigma = sqrt(sigmasq);
  real<lower=0> gamma = sqrt(gammasq);
  vector<lower=0>[p] tau = sqrt(tausq);
  vector<lower=0>[p] betasd = sigma*gamma*tau;
  vector[n] mu;
  mu = beta0 + x*beta;
}

model {
  // priors
  target += jeffreys_lpdf(sigmasq);   
  target += igauss_lpdf(gammasq | ctalphainv, lambda);
  for(j in 1:p) {
    target += igauss_lpdf(tausq[j] | ztthetainv, kappa);
  }
  beta ~ normal(0., betasd);
  
  // likelihood
  y ~ normal(mu, sigma);
}
