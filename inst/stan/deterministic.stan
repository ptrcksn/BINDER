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

 vector<lower=0>[M] psi;
}

transformed parameters {
 vector[N] trans_theta;
 vector<lower=0,upper=1>[N] theta;

 trans_theta = zeta + (X*tau);
 theta = inv_logit(trans_theta);
}

model {
 // Prior:
 zeta ~ normal(mu_zeta,sigma_zeta);
 tau ~ normal(mu_tau,sigma_tau);

 psi ~ cauchy(mu_psi,sigma_psi);

 // Likelihood:
 for(m in 1:M){
  trans_Y[,m] ~ normal(trans_theta,psi[m]);
 }
}
