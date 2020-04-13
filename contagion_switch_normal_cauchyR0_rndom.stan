data {
int T;  /// dimension (time)
int G;  /// nb clusters
int<lower=1> K;  /// latent dimension
real N[G, T];  //data
row_vector<lower=0.0>[K] w;
real alpha0;
real beta0;
real<lower=1>  Di;
}

transformed data {
  real log_unif;
  real NN[G, T];
  log_unif = -log(T);
  for (g in 1:G){
     for (t in 1:T){
         NN[g,t] = 2 * sqrt(N[g,t] +3.0/8.0);
     }
  }
}

parameters{
  real<lower=0.0, upper=1> tau;
  real<lower=0.01> cbar1[G];
  row_vector<lower=0.0>[K] mu0[G];
  row_vector<lower=0.0>[T] r0[G];  
//row_vector<lower=0.0>[G] inv_epsilon;
}

transformed parameters{
    real<lower=0.0> R0[G];
    real<lower=0.0> cbar[G];
    real<lower=0.0>  lambda[G, T];
    vector[T] lp[G];
    //row_vector<lower=0.0>[G] epsilon;
  
    for (g in 1:G){ 
         vector[T + 1] lp_e;
         vector[T + 1] lp_l; 
         cbar[g] = fabs(cbar1[g]);
         R0[g] = cbar[g] * tau * Di;
         //epsilon[g] = 1.0/sqrt(inv_epsilon[g]);
         lambda[g, 1] =  NN[g,1]; /// Lamentable
		 for (t in 2:(K)){    //	
		   lambda[g, t] =  2 * sqrt(3.0/8. + r0[g,t] * dot_product(w[1:(K+1-t)], mu0[g, t:K]) + R0[g] * dot_product(w[(K +2 -t):K], to_vector(N[g, 1:(t-1)]))); 		  
		 }
		 for (t in (K+1):T){    //	
		   lambda[g, t] = 2 * sqrt(3.0/8. + r0[g,t] * dot_product(w , to_vector(N[g, (t- K):(t-1)])));
		 }
		 
        lp_e[1] = 0;
        lp_l[1] = 0;
        for (t in 1:T){
          lp_e[t + 1] = lp_e[t] + normal_lpdf(NN[g, t] | 2.35, 1.0); // early rate Little Poisson
          lp_l[t + 1] = lp_l[t] + normal_lpdf(NN[g, t] | lambda[g,t], 1.0);//epsilon[g]); //  late rate Big Poisson
        }
        lp[g] = rep_vector(log_unif + lp_l[T + 1], T) + head(lp_e, T) - head(lp_l, T);
    }
}

model{
    cbar1[1] ~ normal(1.0, 0.0001);
    for (g in 2:G){ 
        cbar1[g] ~ gamma(alpha0, beta0);
        //inv_epsilon[g] ~ gamma(1.0, 1.0);
        for (t in 1:T){
           r0[g,t] ~ cauchy(R0[G], 1.0);
        }
    }
    
    for (g in 1:G){
    
        for (k in 1:K){ 
            mu0[g, k] ~ gamma(50, 1);
        }  

       target += log_sum_exp(lp[g]);
    } 
   tau ~ beta(1, 39); 
}  
