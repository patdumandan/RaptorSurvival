#final model: beta binomial###########

#use subset of data in Newton database
#fixed effects: mass, type of estimate (direct/indirect survival based on methods used)
#species and study as random effect
#issue: low ESS

#run paralle chains
require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#data#######
surv=read.csv(file.choose(), h=T)
str(surv)
plot(surv$survival.est~surv$Average.mass..kg.)

#standardize data
surv$mass=(surv$Average.mass..kg.-mean(surv$Average.mass..kg.))/(2*sd(surv$Average.mass..kg.))
surv$est.type<-ifelse(surv$death.type=="direct", 0, 1) #apparent (1) /true (0) survival estimate
surv$estimate=(surv$est.type-mean(surv$est.type))/(2*sd(surv$est.type))

#rename data#######

N=length(surv$Species) #no.of rows
surv$spcode=as.integer(surv$Species)
surv$stcode=as.integer(surv$Reference)
surv$famcode=as.integer(surv$family)
Nsp=36 #no.of species
Nst=65 #no of studies
Nfam=6


surv_mod_spstfam=stan(model_code="
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
  vector  [N]death_type;// direct/indirect

                }
 parameters {

  real  alpha;// global intercept
  real alpha_sp[Nsp]; //random intercept per species
  real alpha_st [Nst];// random intercept per study
  real alpha_fam [Nfam];// random intercept per family
  real  beta1; //slope mass
  real  beta2; //slope indirect effect
  real<lower=0> sigma_sp;//errors for random effects
  real<lower=0> sigma_st;//errors for random effects
  real<lower=0> sigma_fam;
  real <lower=0> phi;
  real <lower=0, upper=1> pred_surv[N] ;
              }
   
     
   transformed parameters{
  vector <lower=0, upper=1> [N] surv_mu; //estimated survival 
  vector <lower=0> [N] A;
  vector <lower=0> [N] B;
  
  //model:
  
  for (i in 1:N){
  
  surv_mu[i]= inv_logit(alpha_sp[species[i]]+alpha_st[study[i]]+alpha_fam[family[i]]+beta1*mass[i]+beta2*death_type[i]);
  }
  
  A = surv_mu * phi;
  B = (1 - surv_mu)* phi;// look into this, if phi is not=1, relationship not hold
  
  }
  
  

 model {
  //priors
  
  beta1~ normal (0,1);
  beta2~ normal (0,1);
  sigma_sp ~normal(0,1);
  sigma_st~ normal(0,1);
  sigma_fam~ normal(0,1);
  
  phi ~normal(7,1);// use info. from beta regression of all juv and adult
  
  //model likelihood:
  
  pred_surv ~ beta(A, B); // survival estimate, beta dist.
  y~binomial(n, pred_surv); //no.of survivors drawn from binomial dist; based on sample size and reported survival estimate
 
  alpha_sp~normal(alpha, sigma_sp);
  alpha_st~normal(alpha, sigma_st);
  alpha_fam~normal(alpha, sigma_fam);
   
  }

generated quantities {
  
  real pred_y [N];//predictions on survival

    pred_y = beta_rng(A, B);
   
  }
", data=list(N=length(surv$Species), y=surv$estimated.survived, n=surv$sample.size, mass=surv$mass, death_type=surv$estimate,
             species=surv$spcode,family=surv$famcode, study=surv$stcode,Nfam=6, Nst=65, Nsp=36), chains=4, iter=3000 , warmup=1000) #, control=list(adapt_delta=0.999, max_treedepth=12))

saveRDS(surv_mod_spstfam, file="meta_survival_spstfam.RDS")
pairs(stan_model=surv_mod7, pars=c("alpha", "beta1", "beta2", "lp__", "energy__", "sigma_sp", "sigma_st", "phi"))
post=rstan::extract(surv_mod_spstfam)$pred_y #predicted survival estimate
length(post)
mean(post) 
mean(surv$survival.est) #0.72

#model output visualization
print(surv_mod_spstfam, pars=c("alpha", "beta1", "beta2", "alpha_sp"))

jpeg("predsplot.jpeg", width = 4, height = 4, units = 'in', res = 300)
matplot(surv$Average.mass..kg.,t(post), type="l", col="grey", xlab="average mass (kg)", ylab="survival estimate", ylim=c(0.0,1.0))
points(surv$survival.est~surv$Average.mass..kg., col="black", pch=19)
preds=read.csv(file.choose(), h=T) #estimated survival of other species
points(preds$survival..mod9.~preds$mass..kg., col="white", pch=19)#
dev.off()

print(surv_mod7)

#model output
print(meta_survival_final_v3, pars=c("alpha", "beta1", "beta2", "alpha_sp"))
