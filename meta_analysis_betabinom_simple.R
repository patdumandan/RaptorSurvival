#beta regression  ###########
#no random effects
#mass as predictor
#dropped studies that did not report sample size


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

data {
     int<lower=2> N;          // number of rows
     int<lower=0> y[N];       // survived
     int<lower=0> n[N];       // total 
     real mass[N];            //predictor 1
}
 parameters {

  real <lower=0> alpha;// global intercept
  real <lower=0> beta1; //slope age
  real <lower=0> phi;
  real <lower=0, upper=1> pred_surv;
              }
   
     
  transformed parameters{
  vector <lower=0, upper=1> [N] y_mu; //estimated survival 
  vector <lower=0> [N] A;
  vector <lower=0> [N] B;
  
  for (i in 1:N){
  
  y_mu[i]= inv_logit(alpha+beta1*mass[i]);
  }
  

  A = y_mu * phi;
  B = (1 - y_mu )* phi;// 
  
  }
  
model {
    phi ~normal(7,100);// use infor. from 10 studies 
    pred_surv ~ beta(A, B);
    y~binomial(n, pred_surv);
}
", data=list(N=N, y=y, n=n, Nsp=Nsp, mass=mass), chains=4, iter=500)



