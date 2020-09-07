#final model: beta binomial###########
#non-centered parameterization

#use subset of data in Newton database
#fixed effects: mass, type of estimate (direct/indirect survival based on methods used)
#allowed intercept to vary per species, family and study 
#allowed slopes to vary per family for each trait (mass, diet, foraging strategy)

#run parallel chains
require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#data#
surv=read.csv(file.choose(), h=T)
str(surv)
plot(surv$survival.est~surv$Average.mass..kg.)

#standardize data
surv$mass=(surv$Average.mass..kg.-mean(surv$Average.mass..kg.))/(2*sd(surv$Average.mass..kg.))
surv$forage=(surv$ForStrat.ground-mean(surv$ForStrat.ground))/(2*sd(surv$ForStrat.ground))
surv$diet=(surv$Diet.Inv-mean(surv$Diet.Inv))/(2*sd(surv$Diet.Inv))
surv$est.type<-ifelse(surv$death.type=="direct", 0, 1) #apparent (1) /true (0) survival estimate
surv$estimate=(surv$est.type-mean(surv$est.type))/(2*sd(surv$est.type))
str(surv)

#other input data

N=length(surv$Species) #no.of rows
surv$spcode=as.integer(surv$EnglishName)
surv$stcode=as.integer(surv$Reference)
surv$famcode=as.integer(surv$family)
Nsp=37 #no.of species
Nst=65 #no of studies
Nfam=6 #no of families


surv_mod_trial_spstfam_fin=stan(model_code="
 data{

  int<lower=0> N; // no.of obs
  int <lower=0> y[N];       // survivors
  int <lower=0>  n[N];       // total 
  vector [N] mass;// ave.mass in kg
  vector [N] diet; //invertebrate diet
  vector [N] forage; //ground foraging strategy
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
  real diet_eff;//slope diet
  real for_eff;// slope foraging strat
  real<lower=0> sigma_sp[Nsp];//errors for random effects
  real<lower=0> sigma_st[Nst];//errors for random effects
  real<lower=0> sigma_fam[Nfam];//errors for random effects
  real <lower=0> phi;
  real <lower=0, upper=1> pred_surv[N] ;//survival per observation
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
  
  alpha_sp[j]= alpha+phi*sigma_sp[j];
  }

   for (k in 1:Nst) {
  
  alpha_st[k]= alpha+phi*sigma_st[k];
   }
 
   for (m in 1:Nfam) {
  
 alpha_fam[m]= alpha+phi*sigma_fam[m];
 // mass_fam[m]=mass_eff+phi*sigma_fam[m];
  //diet_fam[m]=diet_eff+phi*sigma_fam[m];
  //for_fam[m]=mass_eff+phi*sigma_fam[m];
  }
 
  
  //model:
  
  for (i in 1:N){
  
  surv_mu[i]= inv_logit(alpha_sp[species[i]]+alpha_st[study[i]]+alpha_fam[family[i]]+
  mass_eff*mass[i]+est_eff*death_type[i]+diet_eff*diet[i]+for_eff*forage[i]);
 // mass_fam[family[i]]+diet_fam[family[i]]+for_fam[family[i]]);
  }
  
  A = surv_mu * phi;
  B = (1 - surv_mu)* phi;
  
  }

 model {
  //priors

  mass_eff~ normal (0.1,1);
  est_eff~ normal (0,1);
  diet_eff~normal(0,1);
  for_eff~normal(0,1);
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
   
  }
", data=list(N=length(surv$Species), y=surv$estimated.survived, n=surv$sample.size,
             diet=surv$diet, forage=surv$forage,
             mass=surv$mass, death_type=surv$estimate,
             species=surv$spcode,family=surv$famcode, study=surv$stcode,
             Nfam=6, Nst=65, Nsp=37), chains=4, iter=3000) # , warmup=1000, control=list(adapt_delta=0.99, max_treedepth=12))

saveRDS(surv_mod_trial_spstfam_fin, file="meta_survival_spstfam_fin.RDS")
post=rstan::extract(surv_mod_trial_spst_rand)$pred_y #predicted survival estimate
print(surv_mod_trial_spstfam_fin, pars=c("alpha", "mass_eff", "est_eff","diet_eff", "for_eff", "alpha_fam"))

#model output visualization
logq=rstan::extract(surv_mod_trial_spstfam_fin)$log_lik
rstan::loo(logq)
matplot(surv$Average.mass..kg.,t(post), type="l", col="grey", xlab="average mass (g)", ylab="survival estimate", ylim=c(0.0,1.0))
points(surv$survival.est~surv$Average.mass..kg., col="black", pch=19)


#sample predictions for species:alpha+mass_eff+est_eff+diet_eff+for_eff+alpha_sp
plogis(-2.7+0.02+-0.06+0.01+-0.11+1.05)

#sample predictions for family:alpha+mass_eff+est_eff+diet_eff+for_eff+alpha_fam
plogis(-2.7+0.02+-0.06+0.01+-0.11+0.28)
