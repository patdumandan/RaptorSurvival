#beta regression  ###########
#fixed effects: mass, type of estimate (direct/indirect survival based on methods used)
#species as random effect

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
Nsp=40
death_type=as.integer(surv$death.type)
str(surv)
study=surv$Reference
Nst=77


surv_mod3=stan(model_code="
  
  data{

  int<lower=0> N; // no.of obs
  real <lower=0, upper=1> survival[N]; //survival estimate
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
  real <lower=0> phi;
              }
   
     
  transformed parameters{
  vector <lower=0, upper=1> [N] surv_mu; //estimated survival 
  vector <lower=0> [N] A;
  vector <lower=0> [N] B;
  
  for (i in 1:N){
  
  surv_mu[i]= inv_logit(alpha+beta1*mass[i]+beta2*death_type[i]+alpha_sp[Nsp]+alpha_st[Nst]);
  }
  

  A = surv_mu * phi;
  B = 1 - surv_mu * phi;
  
  }
  
  model {
  //priors
  
  alpha~ normal (0,2.5);
  beta1~ normal (0,2.5);
  beta2~ normal (0,2.5);
  sigma_sp ~cauchy(0,2.5);
  sigma_st~ cauchy(0,2.5);
  phi~ normal(0,2.5);


  
  for (i in 1:N){
  
  survival[i]~ beta(A, B);

  }

  for(j in 1:Nsp){
           alpha_sp[j]~normal(alpha, sigma_sp);
  }
  
  for (f in 1: Nst){
  
          alpha_st[f]~normal(alpha, sigma_st);
  }
  }
  generated quantities {
  
  real log_lik [N];//predictions
  

    log_lik = beta_rng(A, B);
   
  }", data=list(N=N, survival=survival, mass=mass,death_type=death_type,
                species=species,study=surv$stcode, Nst=77, Nsp=40), chains=2, iter=100)


#model plot########
PS3=readRDS("surv_model_ranspst.RDS")
saveRDS(surv_mod3, file="surv_model_ranspst.RDS")
post3=rstan::extract(surv_model_ranspst)$log_lik
mean(post3)
sp1=post2[,which(surv$spcode==1)] #matrix of one species
matplot(post3,type="l",col="grey",lty=1,xaxt="n", ylab = "survival estimate")
mean_sp3=apply(post3,2,mean)
lines(mean_sp3~surv$Average.mass..kg.,col="red", lwd=2)
sp11=subset(surv, spcode=="1")
points(surv$survival.est, pch=21, cex=1, col="red", bg="red")
plot(mean_sp3~surv$Average.mass..kg.)
plot(surv$survival.est~surv$Average.mass..kg.)
mean(surv$survival.est) #0.765
mean(mean_sp3) #0.637
