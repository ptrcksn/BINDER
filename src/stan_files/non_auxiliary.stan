data {
 int N;
 int M;

 matrix<lower=0,upper=1>[N,M] Y;

 vector<lower=0>[M] mu_psi;
 vector<lower=0>[M] sigma_psi;
}

transformed data{
 matrix[N,M] logit_Y;

 logit_Y = logit(Y);
}

parameters {
 vector[N] logit_theta;

 vector<lower=0>[M] psi;
}

model {
 // Prior:
 psi ~ cauchy(mu_psi,sigma_psi);

 // Likelihood:
 for(m in 1:M){
  logit_Y[,m] ~ normal(logit_theta,psi[m]);
 }
}

generated quantities{
 vector[N] theta;
 
 theta = inv_logit(logit_theta);
}
