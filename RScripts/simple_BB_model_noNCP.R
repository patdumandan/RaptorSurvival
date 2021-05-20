m1=stan(model_code="

data {
  int<lower=0> N ; //number of data points or number of rows in input dataset
  int<lower=0> y[N]; // observed data
  int<lower=0> n[N];
  vector [N] mass;
}

parameters {
  real mu;//intercept
  real<lower = 0> tau;//NCP error
  vector [N] pred_surv;//
  real mass_eff;
}

transformed parameters {
  
  vector [N] surv_mu = mu+mass*mass_eff+pred_surv*tau;
}

model {
  mu ~ normal(0, 1);
  tau ~ normal(0, 1);
  
  pred_surv ~ normal(0, 1);
  
  for(i in 1:N){
    y[i] ~ binomial_logit(n[i], surv_mu[i]);
  }
}

generated quantities {
  vector[N] p = inv_logit(surv_mu);
}",
data=dat_list, chains=2, iter=300)

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

post=extract(m1)$p
matplot(surv$Average.mass..kg.,t(post), type="l",col="grey", xlab="average mass (kg)", ylab="survival estimate")
points(surv$survival.est~surv$Average.mass..kg., col="black", pch=19)

print(m1, pars=c("mass_eff"))
