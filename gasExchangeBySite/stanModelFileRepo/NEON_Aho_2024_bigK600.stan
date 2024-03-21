  
data {
  int<lower = 1> N;
  int<lower = 1> nexpt;
  int<lower = 1> exptID[N];
  vector[N] dist;
  vector[N] logSf6;
  vector[nexpt] Q;
  vector[nexpt] V;
  vector [nexpt] temp;/// one value for each stream temp to convert to K600
  vector[nexpt] w; 
  
}



parameters {
  vector [nexpt] logK600;
  real intercept; //very close to 0
  real<lower = 0> sigma; //  standard deviation
  real a; // real <lower = 0> a;
  real b;
  real <lower = 0> sigma_expt;
  
  
}


transformed parameters{
  
  vector <lower=0> [nexpt] Kd;
  vector <lower=0> [nexpt] Ksf6;
  
  for (j in 1: nexpt){  // make a loop here
  Ksf6[j] = exp(logK600[j]) /  ((600/(3255.3-(217.13*temp[j])+(6.837*temp[j]^2)-(0.08607*temp[j]^3)))^-0.5); //I think we need to divide by z here.. *Q[j]/(V[j]*w[j])
  }
  
  Kd= Ksf6 ./ V;  //Ksf6 (1/day) / V (m/day) = Kd (1/m)
  
  
}

model {
  for (i in 1:N){
    logSf6[i] ~ normal(intercept + -Kd[exptID[i]]*dist[i], sigma); // likelihood
  }
  
  for (j in 1:nexpt){
    logK600[j]~normal( a + b*log(Q[j]) , sigma_expt);
  }
  
  a ~ normal(8.3,10); //0.6
  b ~ normal(-0.3,1); //0.3
  sigma_expt ~ normal (0,2);
  sigma ~ normal (0,0.2);  // added prior on sigma
  intercept ~ normal (0, 0.1); // working with proportions so strong prior on intercept
}

generated quantities {
  vector[N] logSf6_tilde;
  for (n in 1:N){
    logSf6_tilde[n] = normal_rng(intercept + -Kd[exptID[n]]*dist[n], sigma);

}
}



