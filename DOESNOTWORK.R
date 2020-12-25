#remove Andean Condor
surv1=surv%>%
  filter(Species!="Vultur gryphus")
write.csv(surv1, "surv1.csv")
surv1=read.csv("surv1.csv")

surv1$mass=(surv1$Average.mass..kg.-mean(surv1$Average.mass..kg.))/(2*sd(surv1$Average.mass..kg.))
surv1$forage=(surv1$ForStrat.ground-mean(surv1$ForStrat.ground))/(2*sd(surv1$ForStrat.ground))
surv1$diet=(surv1$Diet.Inv-mean(surv1$Diet.Inv))/(2*sd(surv1$Diet.Inv))
surv1$est.type<-ifelse(surv1$death.type=="direct", 0, 1) #apparent (1) /true (0) survival estimate
surv1$spcode=as.integer(surv1$EnglishName)
surv1$stcode=as.integer(surv1$Reference)
surv1$famcode=as.integer(surv1$family)

dat_list=list(
  N=length(surv1$Species), 
  y=surv1$estimated.survived, 
  n=surv1$sample.size,
  diet=surv1$diet, 
  forage=surv1$forage,
  mass=surv1$mass, 
  death_type=surv1$est.type,
  species=surv1$spcode,
  family=surv1$famcode, 
  study=surv1$stcode,
  Nfam=length(unique(surv1$famcode)), 
  Nst=length(unique(surv1$stcode)), 
  Nsp=length(unique(surv1$spcode)))

surv_mod_2=stan(model_code = "

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
  //vector [Nfam] mass_fam; //random slope per family for mass effect
  //vector [Nfam] diet_fam;//random slope per family for diet effect
  //vector [Nfam] for_fam;//random slope per family for foraging effect
  
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
  mass_eff*mass[i]+est_eff*death_type[i]+for_eff*forage[i]+diet_eff*diet[i]);
  }
  
  A = surv_mu * phi;
  B = (1 - surv_mu)* phi;
  
          }

model {
  //priors

  mass_eff~ normal (0.1,1);//mean informed from Newton paper (mass effect)
  est_eff~ normal (0,1);
  diet_eff~normal(0,1);
  for_eff~normal(0,1);
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
   
            }

", data= dat_list, chains= 4, iter=3000, control=list(adapt_delta=0.99))

print(surv_mod_2, pars=c("alpha_fam", "alpha", "mass_eff", "for_eff", "diet_eff"))
saveRDS(surv_mod_2, "survival_all_traits.RDS")
