functions {
  // inverse Gaussian log probability density function
  real igauss_lpdf(real x, real phi, real lambda) {
    return -log(2*pi())/2 + log(lambda)/2 - 3*log(x)/2
           - lambda*square(x - phi)/(2*square(phi)*x);
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
  vector<lower=0>[p] phi;
  real<lower=0> lambdaf;
  real<lower=0> chi;
  real<lower=0> lambdad;
}

parameters {
  vector[p] beta;
  vector<lower=0>[p] gammasq;
  real<lower=0> tausq;
  real<lower=0> sigmasq;
}

transformed parameters {
  vector<lower=0>[p] gamma = sqrt(gammasq);
  real<lower=0> tau = sqrt(tausq);
  real<lower=0> sigma = sqrt(sigmasq);
  vector<lower=0>[p] betasd = sigma*gamma*tau;
  // vector[n] mu;
  // mu = x*beta;
}

model {
  // priors
  target += jeffreys_lpdf(sigmasq);   
  target += igauss_lpdf(tausq | chi, lambdad);
  for(j in 1:p) {
    target += igauss_lpdf(gammasq[j] | phi[j], lambdaf);
  }
  beta ~ normal(0., betasd);
  
  // likelihood
  y ~ normal(x*beta, sigma);
}
