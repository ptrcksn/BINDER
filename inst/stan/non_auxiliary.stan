data {
 int N;
 int M;

 matrix<lower=0,upper=1>[N,M] Y;

 vector<lower=0>[M] mu_psi;
 vector<lower=0>[M] sigma_psi;
}

transformed data{
 matrix[N,M] trans_Y;

 trans_Y = logit(Y);
}

parameters {
 vector[N] trans_theta;

 vector<lower=0>[M] psi;
}

transformed parameters {
 vector<lower=0,upper=1>[N] theta;

 theta = inv_logit(trans_theta);
}

model {
 // Prior:
 psi ~ cauchy(mu_psi,sigma_psi);

 // Likelihood:
 for(m in 1:M){
  trans_Y[,m] ~ normal(trans_theta,psi[m]);
 }
}
