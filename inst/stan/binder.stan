data {
 int N;
 int K;
 int M;

 matrix[N,K] X;
 matrix<lower=0,upper=1>[N,M] Y;


 real mu_zeta;
 real<lower=0> sigma_zeta;
 vector<lower=0>[K] mu_tau;
 vector<lower=0>[K] sigma_tau;

 real<lower=0> mu_phi;
 real<lower=0> sigma_phi;

 vector<lower=0>[M] mu_psi;
 vector<lower=0>[M] sigma_psi;
}

transformed data{
 matrix[N,M] trans_Y;

 trans_Y = logit(Y);
}

parameters {
 real zeta;
 vector<lower=0>[K] tau;

 vector[N] raw_gamma;
 real<lower=0> phi;

 vector<lower=0>[M] psi;
}

transformed parameters {
 vector[N] gamma;
 vector[N] trans_theta;
 vector<lower=0,upper=1>[N] theta;

 gamma = zeta + (X*tau);
 trans_theta = (raw_gamma + gamma) * phi;
 theta = inv_logit(trans_theta);
}

model {
 // Prior:
 zeta ~ normal(mu_zeta,sigma_zeta);
 tau ~ normal(mu_tau,sigma_tau);

 raw_gamma ~ normal(0,1);
 phi ~ normal(mu_phi,sigma_phi);

 psi ~ normal(mu_psi,sigma_psi);

 // Likelihood:
 for(m in 1:M){
  trans_Y[,m] ~ normal(trans_theta,psi[m]);
 }
}
