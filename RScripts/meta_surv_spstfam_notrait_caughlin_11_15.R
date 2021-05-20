#next steps:
#have two separate versions of code for singletons vs those that have multiple individuals
#tytonidae, pandionidae, 	Cathartidae: all have one species represented

#actually, pull the predictions out of stan. that should let you do everything, I think

recodeFACTOR<-function(FACTOR) {
  number<-as.numeric(as.factor(as.numeric(FACTOR)))
  should<-c(1:max(number))
  if(length(which(should %in% number==F))>0) {
    return("BROKE BROKE BROKE")
  }
  else {return(number)}
}

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
setwd("~/Dropbox/MISC academic projects/raptor meta")
############

surv=read.csv("raptor_survival_cleaned.csv")
str(surv)
plot(surv$survival.est~surv$Average.mass..kg.)

#standardize data
surv$mass=(surv$Average.mass..kg.-mean(surv$Average.mass..kg.))/(2*sd(surv$Average.mass..kg.))
surv$forage=(surv$ForStrat.ground-mean(surv$ForStrat.ground))/(2*sd(surv$ForStrat.ground))
surv$diet=(surv$Diet.Inv-mean(surv$Diet.Inv))/(2*sd(surv$Diet.Inv))
surv$est.type<-ifelse(surv$death.type=="direct", 0, 1) #apparent (1) /true (0) survival estimate
str(surv)

surv$EnglishName<-as.factor(trimws(surv$EnglishName))
surv$Species<-as.factor(trimws(surv$Species))
surv$family<-as.factor(trimws(surv$family))
surv$Reference<-as.factor(trimws(surv$Reference))

surv<-surv[which(surv$EnglishName!="Andean Condor"),]


#other input data

dat_list_train<-vector("list",length=length(unique(surv$Species)))

dat_list_test<-vector("list",length=length(unique(surv$Species)))

for(i in 1:length(unique(surv$Species))){

surv_sub<-subset(surv,surv$Species!=unique(surv$Species)[i])
surv_out<-subset(surv,surv$Species==unique(surv$Species)[i])

key=data.frame(surv_sub$family, recodeFACTOR(droplevels(surv_sub$family)))

keyd=which(surv_sub$family %in% surv_out$family)[1]

surv_out$famcode=rep(keyd,times=length(surv_out$family))

dat_list_train[[i]]=list(
#training data
N=length(surv_sub$Species), y=surv_sub$estimated.survived, n=surv_sub$sample.size,
             diet=surv_sub$diet, forage=surv_sub$forage,
             mass=surv_sub$mass, death_type=surv_sub$est.type,
             species=recodeFACTOR(droplevels(surv_sub$Species)),
              family=recodeFACTOR(droplevels(surv_sub$family)),
                
study=recodeFACTOR(droplevels(surv_sub$Reference)),
             Nfam=length(unique(surv_sub$family)), 
Nst=length(unique(surv_sub$Reference)), Nsp=length(unique(surv_sub$Species)),
             
#test data             
N_test=length(surv_out$Species), n_test=surv_out$sample.size,
             diet_test=surv_out$diet, forage_test=surv_out$forage,
             mass_test=surv_out$mass, death_type_test=surv_out$est.type,
             family_test=surv_out$famcode, 
             Nfam_test=length(unique(surv_out$family)),Nst_test=length(unique(surv_out$Reference)), 
Nsp_test=length(unique(surv_out$Species))
           
                   )


dat_list_test[[i]]=list(
  #test data
  y_real=surv_out$estimated.survived,
  n_test=surv_out$sample.size
)

}


surv_mod_notrait_spstfam_fin=stan(file="bird_leave_one_out.stan", 
data=dat_list_train[[2]], chains=4, iter=3000,
pars=c("surv_mu_test","phi","pred_y_test","y_test","mass_eff","est_eff")) # , warmup=1000, control=list(adapt_delta=0.99, max_treedepth=12))
#

#two problems here:
#1) when family is just one it breaks
#2)   Stan does not support NA (in family_test) in data
table(surv$family)


allpreds=vector("list",length(dat_list_test))

#when it is just one it breaks
op=34
dat_list_train[[op]]$family_test #why?

for(op in 1:length(dat_list_test)){
  ##
  ##
  ##
  ##
  ##put in separate indicator here for models with just ONE row in test data
  runit=stan(file="bird_leave_one_out.stan", 
     data=dat_list_train[[op]], chains=4, iter=3000,
    pars=c("surv_mu_test","phi","pred_y_test","y_test","mass_eff","est_eff")) # , warmup=1000, control=list(adapt_delta=0.99, max_treedepth=12))
  #
  rundat=as.data.frame(runit)
  
  
  allpreds[[op]]<-rundat[,grep("y_test",colnames(rundat))]
  
  print(op)
  
}

save(allpreds,file="predictions_first_demo.Rdata")
#
runbad=stan(file="bird_leave_one_out.stan", 
           data=dat_list_train[[4]], chains=4, iter=3000,
           pars=c("surv_mu_test","phi","pred_y_test","y_test","mass_eff","est_eff")) # , warmup=1000, control=list(adapt_delta=0.99, max_treedepth=12))
#


saveRDS(surv_mod_notrait_spstfam_fin, file="meta_survival_spstfam_notrait.RDS")
post=rstan::extract(surv_mod_notrait_spstfam_fin)$pred_y #predicted survival estimate
print(surv_mod_notrait_spstfam_fin, pars=c("alpha", "mass_eff", "alpha_fam"))

#model output visualization
logq=rstan::extract(surv_mod_notrait_spstfam_fin)$log_lik
rstan::loo(logq)
matplot(surv$Average.mass..kg.,t(post), type="l", col="grey", xlab="average mass (g)", ylab="survival estimate", ylim=c(0.0,1.0))
points(surv$survival.est~surv$Average.mass..kg., col="black", pch=19)


#sample predictions for species:alpha+mass_eff+est_eff+diet_eff+for_eff+alpha_sp
plogis(-2.7+0.02+-0.06+0.01+-0.11+1.05)

#sample predictions for family:alpha+mass_eff+est_eff+diet_eff+for_eff+alpha_fam
plogis(-2.7+0.02+-0.06+0.01+-0.11+0.28)
