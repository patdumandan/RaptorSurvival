new_mod=stan(model_code = "
             data{

  int<lower=0> N; // no.of obs
  int <lower=0> y[N];       // survivors
  int <lower=0>  n[N];       // total 
  vector [N] mass;// ave.mass in kg
  int species[N]; //ID of each species
  int family [N]; //ID of family
  int study [N]; //ID of study
  int Nsp; //no.of species
  int Nst; //no.of studies
  int Nfam;// no. of families
  vector[N] death_type;// direct/indirect
 }
                
parameters {

  real alpha;// global intercept
  real mass_eff; //slope mass
  real est_eff; //slope indirect effect
  real<lower=0> sigma_sp[Nsp];//errors for random effects
  real<lower=0> sigma_st[Nst];//errors for random effects
  real<lower=0> sigma_fam[Nfam];//errors for random effects
  real <lower=0> phi;//beta error term
  real <lower=0> sp_non;//non-centered error term for species
  real <lower=0> st_non;//non-centered error term for study
  real <lower=0> fam_non;//non-centered error term for family
  real <lower=0, upper=1> pred_surv[N] ;//survival per observation
}
   
transformed parameters{
  
  vector <lower=0, upper=1> [N] surv_mu; //mean estimated survival 
  vector <lower=0> [N] A;
  vector <lower=0> [N] B;
  vector [Nsp] alpha_sp; //random intercept per species
  vector [Nst] alpha_st;// random intercept per study
  vector [Nfam] alpha_fam;// random intercept per family

  for (j in 1:Nsp) {
  
  alpha_sp[j]= sp_non*sigma_sp[j];
  }

   for (k in 1:Nst) {
  
  alpha_st[k]= st_non*sigma_st[k];
   }
 
   for (m in 1:Nfam) {
  
 alpha_fam[m]= fam_non*sigma_fam[m];
 
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
  sp_non~ normal(0,1);
  st_non~ normal(0,1);
  fam_non~ normal(0,1);
  sigma_sp~ normal(0,1);
  sigma_st~ normal(0,1);
  sigma_fam~ normal(0,1);
  
  phi ~normal(7,1);// use info. from beta regression of all juv and adult
  
  //model likelihood:
  
  pred_surv ~ beta(A, B); // survival estimate, beta dist.
  y~binomial(n, pred_surv); //no.of survivors drawn from binomial dist; based on sample size and reported survival estimate
 
  }

generated quantities {
  
  real pred_y [N];//predictions on survival
  real log_lik [N];// for looic calculations
  
    pred_y = beta_rng(A, B);
    
    for (x in 1:N){
    log_lik[x]= beta_lpdf(pred_surv[x]| A[x], B[x]);}
   
  }",data=data_list,chains=4, iter=3000, control=list(adapt_delta=0.99))
  

data_list=list(N=length(surv$Species), y=surv$estimated.survived, n=surv$sample.size,
               diet=surv$diet, forage=surv$forage,
               mass=surv$mass, death_type=surv$estimate,
               species=surv$spcode,family=surv$famcode, study=surv$stcode,
               Nfam=length(unique(surv$family)), Nst=length(unique(surv$Reference)), Nsp=length(unique(surv$EnglishName))) 

print(new_mod, pars=c("alpha", "alpha_sp"))
saveRDS(new_mod, "new_model.RDS")
