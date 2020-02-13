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
  int<lower=1> D;
  int p[D];
  int<lower=1> sump;
  int<lower=1> n;
  int<lower=1> G;
  int<lower=1> H;
  matrix[n, sump] x;
  matrix[n, D] y; 
  matrix[sump, G] C;
  matrix[D, H] Z;
  real<lower=0> nuf;
  real<lower=0> nud;
  real<lower=0> kappaf;
  real<lower=0> kappad;
  real<lower=0> xif;
  real<lower=0> xid;
}

parameters {
  vector[sump] beta;
  vector<lower=0>[sump] gammasq;
  vector<lower=0>[D] tausq;
  vector<lower=0>[D] sigmasq;
  vector[G] alphaf;
  vector[H] alphad;
  real<lower=0> lambdaf;
  real<lower=0> lambdad;
}

transformed parameters {
  vector<lower=0>[sump] gamma = sqrt(gammasq);
  vector<lower=0>[D] tau = sqrt(tausq);
  vector<lower=0>[D] sigma = sqrt(sigmasq);
  vector<lower=0>[sump] phi = C*alphaf;
  vector<lower=0>[D] chi = Z*alphad;
}

model {
  int cp;
  int pos;
  pos = 1;
  for(d in 1:D) {
    target += jeffreys_lpdf(sigmasq[d]);   
    target += igauss_lpdf(tausq[d] | chi[d], lambdad);
    
    cp = p[d];
    for(j in 1:cp) {
      target += igauss_lpdf(gammasq[pos + j - 1] | phi[pos + j - 1], lambdaf);
      beta[pos + j - 1] ~ normal(0., gamma[pos + j - 1]*tau[d]*sigma[d]);
    }
    y[, d] ~ normal(block(x, 1, pos, n, p[d])*segment(beta, pos, p[d]), sigma[d]);
    pos = pos + p[d];
  }
  alphaf ~ normal(0., nuf/lambdaf);
  alphad ~ normal(0., nud/lambdad);
  lambdaf ~ gamma(kappaf, lambdaf);
  lambdad ~ gamma(kappad, lambdad);
}
