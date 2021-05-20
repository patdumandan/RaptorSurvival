//severe underprediction
data{

  int<lower=0> N; // no.of obs
  int <lower=0> y[N];       // survivors
  int <lower=0>  n[N];       // total 
  vector [N] mass;// ave.mass in kg
//  vector [N] diet; //invertebrate diet
//  vector [N] forage; //ground foraging strategy
  int species[N]; //ID of each species
  int family [N]; //ID of family
  int study [N]; //ID of study
  int Nsp; //no.of species
  int Nst; //no.of studies
  int Nfam;// no. of families
  vector[N] death_type;// direct/indirect

//test data

  int<lower=0> N_test; // no.of obs
  int <lower=0>  n_test;       // total 
  real mass_test;// ave.mass in kg
//  vector [N_test] diet_test; //invertebrate diet
//  vector [N_test] forage_test; //ground foraging strategy
  int<lower=0> family_test; //ID of family
  int Nsp_test; //no.of species
  int Nst_test; //no.of studies
  int Nfam_test;// no. of families
  int<lower=0> death_type_test;// direct/indirect

 }
                
 parameters {

  real alpha;// global intercept
  real mass_eff; //slope mass
  real est_eff; //slope indirect effect
//  real diet_eff;//slope diet
//  real for_eff;// slope foraging strat
  real<lower=0> sigma_sp;//errors for random effects
  real<lower=0> sigma_st;//errors for random effects
  real<lower=0> sigma_fam;//errors for random effects
  real <lower=0> phi;
  real <lower=0> sp_non;//non-centered error term for species
  real <lower=0> st_non;//non-centered error term for study
  real <lower=0> fam_non;//non-centered error term for family
  real <lower=0, upper=1> pred_surv[N] ;//survival per observation
  real <lower=0, upper=1> pred_surv_test[N_test];
              }
   
     
   transformed parameters{
  vector <lower=0, upper=1> [N] surv_mu; //mean estimated survival 
  vector <lower=0> [N] A;
  vector <lower=0> [N] B;
  vector [Nsp] alpha_sp; //random intercept per species
  vector [Nst] alpha_st;// random intercept per study
  vector [Nfam] alpha_fam;// random intercept per family
  //vector [Nfam] mass_fam; //random slope per family for mass effect
  //vector [Nfam] diet_fam;//random slope per family for diet effect
  //vector [Nfam] for_fam;//random slope per family for foraging effect
  
  for (j in 1:Nsp) {
  
  alpha_sp[j]= sp_non*sigma_sp;
  }
   for (k in 1:Nst) {
  
  alpha_st[k]= st_non*sigma_st;
   }
 
   for (m in 1:Nfam) {
  
 alpha_fam[m]= fam_non*sigma_fam;
 
  }
  
  //model:
  
  for (i in 1:N){
  
  surv_mu[i]= inv_logit(alpha+alpha_sp[species[i]]+alpha_st[study[i]]+alpha_fam[family[i]]+
  mass_eff*mass[i]+est_eff*death_type[i]);
  }
  
  A = surv_mu * phi;
  B = (1 - surv_mu)* phi;
  
  }

 model {
  //priors

  mass_eff~ normal (0.1,1);
  est_eff~ normal (0,1);
  //diet_eff~normal(0,1);
  //for_eff~normal(0,1);
  sp_non~ normal(0,1);
  st_non~ normal(0,1);
  fam_non~ normal(0,1);
  sigma_sp~ normal(0,1);
  sigma_st~ normal(0,1);
  sigma_fam~ normal(0,1);
  alpha~normal(1,0.5);
  
  phi ~normal(7,1);// use info. from beta regression of all juv and adult
  
  //model likelihood:
  
  pred_surv ~ beta(A, B); // survival estimate, beta dist.
  y~binomial(n, pred_surv); //no.of survivors drawn from binomial dist; based on sample size and reported survival estimate
 
  }

generated quantities {
  
    real <lower=0, upper=1> surv_mu_test;
   real <lower=0>  A_test;
  real <lower=0>  B_test;
   real pred_y_test;//predictions on survival
 int<lower=0> y_test; //predictions on data level
 

  
   
  //surv_mu[i]= inv_logit(alpha+alpha_sp[species[i]]+alpha_st[study[i]]+alpha_fam[family[i]]+
  //mass_eff*mass[i]+est_eff*death_type[i]);
  
  surv_mu_test= inv_logit(alpha+alpha_fam[family_test]+
  mass_eff*mass_test+est_eff*death_type_test);
  //+diet_eff*diet_test[ii]+for_eff*forage_test[ii]);
 // mass_fam[family_test[ii]]+diet_fam[family_test[ii]]+for_fam[family_test[ii]]);
  
  A_test = surv_mu_test * phi;
  B_test = (1 - surv_mu_test)* phi;
  
    pred_y_test = beta_rng(A_test, B_test);
    y_test= binomial_rng(n_test, pred_y_test); 
}