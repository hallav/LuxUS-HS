data {

  int<lower=1> n_cytosines;
  int<lower=1> n_replicates;
  int<lower=1> n_predictors;

  real<lower=0,upper=1> bsEff[n_replicates * n_cytosines];
  real<lower=0,upper=1> bsBEff[n_replicates * n_cytosines];
  real<lower=0,upper=1> seqErr[n_replicates * n_cytosines];

  int<lower=0> bsC[n_replicates * n_cytosines];
  int<lower=0> bsTot[n_replicates * n_cytosines];

  matrix[n_cytosines * n_replicates, n_cytosines * n_predictors] X;

  int Z_R[n_replicates * n_cytosines];

  real<lower=0> sigmaB2;

  real<lower=0> alpha;
  real<lower=0> beta;

  real<lower=0> alphal;
  real<lower=0> betal;

  real<lower=0> alphaR;
  real<lower=0> betaR;

  real<lower=0> coordinates[n_cytosines * n_predictors];

}
parameters {

  vector[n_cytosines * n_predictors] B;
  real<lower=0> sigmaE2;
  real<lower=0> sigmaR2;
  vector[n_cytosines * n_replicates] Y;
  real<lower=0> l;
  vector[n_replicates] ranint_rep;
  vector<lower=0>[n_cytosines] z;
  real<lower=0> r1_global;
  real<lower=0> r2_global;
  real<lower=0> r1_local[n_cytosines];
  real<lower=0> r2_local[n_cytosines];

}
transformed parameters {

  vector<lower=0,upper=1>[n_cytosines * n_replicates] theta;
  cov_matrix[n_cytosines * n_predictors] SigmaB; 
  vector<lower=0,upper=1>[n_cytosines * n_predictors] d; 
  real<lower=0> tau;
  real<lower=0> lambda[n_cytosines]; 
  real d_tilde[n_cytosines];
  vector[n_cytosines * n_replicates] Y_hat;

  tau = r1_global * sqrt(r2_global);
  
  for (i in 1:n_cytosines){
    lambda[i] = r1_local[i] * sqrt(r2_local[i]);
    d_tilde[i] = z[i] * lambda[i] * tau;
  } 
  
  for (i in 1:n_cytosines){ 
    for (j in 1:n_predictors){
      d[n_predictors * (i - 1) + j] = 1 - 1 / pow(1 + 10 * exp(-5 * d_tilde[i]), 1/0.5);
    }
  }


  for (i in 1:n_predictors * n_cytosines){
    for (j in 1:n_predictors * n_cytosines){
      if (i != j && fmod(i,n_predictors)==fmod(j,n_predictors)) {
        SigmaB[i,j] = exp(-(fabs(coordinates[j]-coordinates[i])) / pow(l,2)) * d[i] * d[j];
      }else{
        SigmaB[i,j] = 0;
      }
    }
  }

  for (i in 1:(n_predictors * n_cytosines)){
    SigmaB[i,i] = 1;
  }


  for (i in 1:n_cytosines * n_replicates){
    Y_hat[i] = dot_product(X[i,:], B) + ranint_rep[Z_R[i]];
  }

  for (i in 1:(n_cytosines * n_replicates)){
    theta[i] = inv_logit(Y[i]);
  }

}
model {

  l ~ gamma(alphal,betal);

  for (i in 1:n_cytosines){
    z[i] ~ normal(0,1) T[0,];
  }
  r1_global ~ normal(0,1);
  r2_global ~ inv_gamma(0.5,0.5);

  for (i in 1:n_cytosines){
    r1_local[i] ~ normal(0,1);
    r2_local[i] ~ inv_gamma(0.5,0.5);
  }

  sigmaE2 ~ gamma(alpha,beta);
  sigmaR2 ~ gamma(alphaR,betaR);

  B ~ multi_normal(rep_vector(0,n_predictors * n_cytosines), sigmaB2 * SigmaB);

  ranint_rep ~ multi_normal(rep_vector(0,n_replicates), diag_matrix(rep_vector(sigmaR2,n_replicates)));

  Y ~ multi_normal(Y_hat, diag_matrix(rep_vector(sigmaE2, n_cytosines * n_replicates)));

  for (i in 1:n_cytosines * n_replicates){
    bsC[i] ~ binomial(bsTot[i],
            theta[i] * ((1.0 - seqErr[i]) * (1.0 - bsEff[i]) + seqErr[i] * bsEff[i]) +
            (1-theta[i]) * ((1.0 - bsBEff[i]) * (1.0 - seqErr[i]) + seqErr[i] * bsBEff[i]));
  }

}
