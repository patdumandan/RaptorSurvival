library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

surv=read.csv("https://raw.githubusercontent.com/patdumandan/RaptorSurvival/master/surv1.csv")

surv$mass=(surv$Average.mass..kg.-mean(surv$Average.mass..kg.))/(2*sd(surv$Average.mass..kg.))
surv$estimate<-ifelse(surv$death_type=="direct", 0, 1) #apparent (1) /true (0) survival estimate
surv$spcode=as.integer(surv$Species)
surv$famcode=as.integer(surv$family)
surv$stcode=as.integer(surv$Reference)

surv<-surv[which(surv$EnglishName!="Andean Condor"),] 
str(surv)

dat_list=list( N=length(surv$Reference),
               y=surv$estimated.survived,    
               n=surv$sample.size,     
               mass=surv$mass,
               species=surv$spcode,
               family=surv$famcode,
               study=surv$stcode,
               Nsp=length(unique(surv$spcode)),
               Nst=length(unique(surv$stcode)),
               Nfam=length(unique(surv$famcode)),
               death_type=surv$estimate)

surv_mod=stan(model_code="

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
  real mass_eff; //slope for mass
  real est_eff; //slope for estimate type

  real<lower=0> sigma; //for NCP
  real <lower=0> sigma_sp[Nsp];//errors for random effects
  real<lower=0> sigma_st[Nst];//errors for random effects
  real<lower=0> sigma_fam[Nfam];//errors for random effects
 
  real <lower=0> phi;// variance of the likelihood
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
  
  alpha_sp[j]= alpha+sigma*sigma_sp[j];
  }

   for (k in 1:Nst) {
  
  alpha_st[k]= alpha+sigma*sigma_st[k];
   }
 
   for (m in 1:Nfam) {
  
 alpha_fam[m]= alpha+sigma*sigma_fam[m];
  }
 
  
  //model:
  
  for (i in 1:N){
  
  surv_mu[i]= inv_logit(alpha_sp[species[i]]+alpha_st[study[i]]+alpha_fam[family[i]]+
              mass_eff*mass[i]+est_eff*death_type[i]);

  }
  
  A = surv_mu * phi;
  B = (1 - surv_mu)* phi;
  
  }

 model {
  //priors

  mass_eff~ normal (0.1,1);
  est_eff~ normal (0,1);
 
  sigma~normal(0,1);
  sigma_sp~ normal(0,1);
  sigma_st~ normal(0,1);
  sigma_fam~ normal(0,1);
  
  phi ~normal(7,1);// use info. from beta regression of all juv and adult
  
  //model likelihood:
  
  pred_surv ~ beta(A, B); // survival estimate, beta dist.
  y~binomial(n, pred_surv); //no.of survivors drawn from binomial dist; based on sample size and reported survival estimate
 
  }", data=dat_list, chains=4, iter=3000)

saveRDS(surv_mod, file="bbmod_ncp.RDS")
post=extract(surv_mod)$pred_surv
matplot(surv$Average.mass..kg.,t(post), type="l",col="grey", xlab="average mass (kg)", ylab="survival estimate")
points(surv$survival.est~surv$Average.mass..kg., col="black", pch=19)

print(surv_mod, pars=c("alpha", "mass_eff", "alpha_sp", "alpha_fam", "alpha_st"))
