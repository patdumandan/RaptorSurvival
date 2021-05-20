#beta binomial###########
#fixed effects: mass, type of estimate (direct/indirect survival based on methods used)
#species and study as random effect


#run paralle chains
mc.cores = parallel::detectCores()

#data#######
surv=read.csv(file.choose(), h=T)
str(surv)
plot(surv$survival.est~surv$Average.mass..kg.)
survi=surv$survival.est
surv$mass=surv$Average.mass..kg.
N=length(surv$Species)
survival=surv$survival.est
mass=surv$Average.mass..kg.
species=surv$spcode
surv$spcode=as.integer(surv$Species)
surv$stcode=as.integer(surv$Reference)
Nsp=36
death_type=as.integer(surv$death.type)
str(surv)
study=surv$Reference
Nst=65
y=surv$total.survived

surv_mod5=stan(model_code="

 data{

  int<lower=0> N; // no.of obs
  int<lower=0> y[N];       // survived
  int<lower=0> n[N];       // total 
  real <lower=0> mass[N];// ave.mass in kg
  int species[N]; //ID of each species
  int study [N]; //ID of study
  int Nsp; //no.of species
  int Nst; //no.of studies
  int death_type[N];// direct/indirect

                }
 parameters {

  real <lower=0> alpha;// global intercept
  real <lower=0> alpha_sp[Nsp]; //random intercept per species
  real <lower=0> alpha_st [Nst];// random intercept per study
  real <lower=0> beta1; //slope age
  real <lower=0> beta2; //slope indirect effect
  real<lower=0> sigma_sp;//errors for random effects
   real<lower=0> sigma_st;//errors for random effects
  real <lower=0> phi;
  real <lower=0, upper=1> pred_surv;
              }
   
     
   transformed parameters{
  vector <lower=0, upper=1> [N] surv_mu; //estimated survival 
  vector <lower=0> [N] A;
  vector <lower=0> [N] B;
  
  for (i in 1:N){
  
  surv_mu[i]= inv_logit(alpha+beta1*mass[i]+beta2*death_type[i]+alpha_sp[Nsp]+alpha_st[Nst]);
  }
  
  A = surv_mu * phi;
  B = (1 - surv_mu )* phi;// look into this, if phi is not=1, relationship not hold
  
  }
  
  

 model {
  //priors
  
  alpha~ normal (0,1);
  beta1~ normal (0,1);
  beta2~ normal (0,1);
  sigma_sp ~normal(0,1);
  sigma_st~ normal(0,1);
  phi ~normal(7,100);// use infor. from 10 studies 
  
  
  pred_surv ~ beta(A, B);
  y~binomial(n, pred_surv);

for(j in 1:Nsp){
           alpha_sp[j]~normal(0, sigma_sp);
  }
  
  for (f in 1: Nst){
  
          alpha_st[f]~normal(0, sigma_st);
  }
  }

generated quantities {
  
  real log_lik [N];//predictions
  

    log_lik = beta_binomial_rng(n, A, B);
   
  }
", data=list(N=N, y=y, n=n, mass=mass, death_type=death_type,
             species=species,study=surv$stcode, Nst=65, Nsp=36), chains=2, iter=2000, control = list(max_treedepth = 12))


#model diagnostics
saveRDS(surv_mod5, file="surv_betabinom_raneff.RDS")
post5=rstan::extract(surv_mod5)$log_lik #predicted number of survivors
mean(post5) #205.4579
mean(surv$total.survived) #205.3369
post51=rstan::extract(surv_mod5)$pred_surv #predicted survival est
matplot(t(post51), type="l", ylim=c(0.4,0.9))
points(surv$survival.est, col="red", lwd=2)
mean(post51) #0.7725
mean(surv$survival.est) #0.7479
