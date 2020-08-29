#final model: beta binomial###########

#use subset of data in Newton database
#fixed effects: mass, type of estimate (direct/indirect survival based on methods used)
#species and study as random effect
#issue: low ESS

#run paralle chains
mc.cores = parallel::detectCores()

#data#######
surv=read.csv(file.choose(), h=T)
str(surv)
plot(surv$survival.est~surv$Average.mass..kg.)

#rename data#######

N=length(surv$Species) #no.of rows
surv$spcode=as.integer(surv$Species)
surv$stcode=as.integer(surv$Reference)

Nsp=36 #no.of species
Nst=65 #no of studies

surv$est.type<-ifelse(surv$death.type=="direct", 0, 1) #apparent (1) /true (0) survival estimate


surv_mod7=stan(model_code="
 data{

  int<lower=0> N; // no.of obs
  int<lower=0> y[N];       // survivors
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
  real <lower=0> beta1; //slope mass
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
  
  //model:
  
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
 
  phi ~normal(7,1);// use info. from beta regression of all juv and adult
  
  //model likelihood:
  
  pred_surv ~ beta(A, B); // survival estimate, beta dist.
  y~binomial(n, pred_surv); //no.of survivors drawn from binomial dist; based on sample size and reported survival estimate

  for(j in 1:Nsp){
           alpha_sp[j]~normal(0, sigma_sp);
  }
  
  for (f in 1: Nst){
  
          alpha_st[f]~normal(0, sigma_st);
  }
  
   
  }

generated quantities {
  
  real log_y [N];//predictions on survival

    log_y = beta_rng(A, B);
   
  }
", data=list(N=N, y=surv$estimated.survived, n=surv$sample.size, mass=surv$Average.mass..kg., death_type=surv$est.type,
             species=surv$spcode,study=surv$stcode, Nst=65, Nsp=36), chains=4, iter=5000, warmup=1000, control=list(adapt_delta=0.99, max_treedepth=12))

saveRDS(surv_mod7, file="meta_survival_final.RDS")
post=rstan::extract(surv_mod7)$log_y #predicted survival estimate
length(post)
mean(post) #0.73 
mean(surv$survival.est) #0.72

#model output visualization
summary(post7)

jpeg("predsplot.jpeg", width = 4, height = 4, units = 'in', res = 300)
matplot(surv$Average.mass..kg.,t(post), type="l", col="grey", xlab="average mass(kg)", ylab="survival estimate", ylim=c(0.1,1.0))
points(surv$survival.est~surv$Average.mass..kg., col="black", pch=19)
preds=read.csv(file.choose(), h=T) #estimated survival of other species
points(preds$mean.estimated.survival.inv_logit.~preds$mass..kg., col="white", pch=19)#
dev.off()
print(surv_mod7)

#model output
print(surv_mod7, pars=c("alpha", "beta1", "beta2"))