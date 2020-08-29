surv_mod9=stan(model_code="

data{
  
  int N; //no.of obs
  int Nsp; //no of species
 // int Nst; //no of studies
  int <lower=1, upper=Nsp> species[N];//species code
//  int <lower=1, upper=Nst> study[N];//study code
  real  mass[N]; //mass kg
  int  y[N]; //no of survivors
  int  n[N]; // sample size
  int death_type [N];// apparent/true estimate
}

parameters{
  vector <lower=0> [Nsp] alpha1; //per species intercept
 // vector <lower=0>[Nst] alpha2; // per study intercept
  real mass_eff;
  real <lower=0> alpha1_mu;
  //real alpha2_mu;
  real <lower=0>sigma_sp;
//  real <lower=0>sigma_st;
  real <lower=0> phi;//variance shape params
  real <lower=0, upper=1> pred_surv;//modelled survival
}

transformed parameters{
  
  vector<lower=0, upper=1> [N] surv_mu;
  vector <lower=0> [N] A;
  vector <lower=0> [N] B;
 
   for (i in 1:N){
    surv_mu[i]=inv_logit(alpha1[species[i]]+(mass[i]*mass_eff));
  
    A = surv_mu * phi;
    B = (1 - surv_mu )* phi;// look into this, if phi is not=1, relationship not hold
   }
}
  model{
    
    sigma_sp~normal(0,1);
    alpha1_mu~normal(0,1);
   // sigma_sp~normal(0,1);
    alpha1~normal(alpha1_mu, sigma_sp);
   // alpha2~normal(alpha2_mu, sigma_st);
    mass_eff~normal(0,1);
    phi~normal(7,1);

    pred_surv ~ beta(A, B); // survival estimate, beta dist.
    y~ binomial(n, pred_surv); //no.of survivors drawn from binomial dist; based on sample size and reported survival estimate
    
  }
  generated quantities {
    
    real y_pred [N];//predictions on survival
    
    y_pred = beta_rng(A, B);
    
  }",

data= list(N= length(surv$Reference), Nsp=36, Nst=65, species=surv$spcode, study=surv$stcode,
mass=surv$Average.mass..kg., y=surv$estimated.survived, n=surv$sample.size, death_type=surv$est.type),
chains=2, iter=200)

