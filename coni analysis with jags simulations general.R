### model calculating the annual indices of abundance from
### nightjar program and the BBS, on paired routes

simparam = data.frame(simparam = c("base","10pdecrease","stable","20yr"),
                        nyearspred = c(11,11,11,21),
                      B.true = c(log(((0.7)^(1/10))),
                                 log(((0.9)^(1/10))),
                                 0,
                                 log(((0.7)^(1/10)))),
                      null1 = c(0,
                                0,
                                log(((0.7)^(1/10))),
                                0),
                      null2 = c(0,
                                log(((0.7)^(1/10))),
                                log(((0.9)^(1/10))),
                                0))

#simparam = simparam[c(1:3),] #removes the 20-year simulation, which showed almost 100% power for all comparisons



library(jagsUI)
library(loo)
library(lme4)

library(reshape2)
library(progress)


for(sj in 1:nrow(simparam)){
set.seed(2019)
  #setwd(paste0("/nightjarfinal"))
  
    sim = simparam[sj,"simparam"]

# data-preparation --------------------------------------------------------


load("full bbs routelist.RData")


load("CONI BBS data.RData")

birds = birds[which(birds$runyear > 2005 & birds$countrynum == 124),]
bbsrtsw4200 = unique(birds[,c("BCR","rt.uni")])



bbsrtsbybcr = table(bbsrtsw4200$BCR)



nj1 = read.csv("targeted_survey_data.csv",
              stringsAsFactors = F)

nj = unique(nj1[,c("region",
                   "BCR",
                   "route",
                   "year")])

bcrs = c(4,6,11)

for(i in 1:nrow(nj)){
  nj[i,"count"] = sum(nj1[which(nj1$region == nj[i,"region"] &
                                  nj1$BCR == nj[i,"BCR"] &
                                  nj1$route == nj[i,"route"] &
                                  nj1$year == nj[i,"year"]),"coni"],na.rm = F)
}

nj = nj[which(nj$BCR %in% bcrs),]

nj$route = paste(nj$BCR,nj$route,sep = "_")

fyearnj = min(nj$year) #first year of the cns survey used as the centering parameter for the slope

bbsalt1 = read.csv("BBS_data_matching.csv")


bbsalt = unique(bbsalt1[,c("region",
                   "BCR",
                   "route",
                   "year")])

for(i in 1:nrow(bbsalt)){
  bbsalt[i,"count"] = sum(bbsalt1[which(bbsalt1$region == bbsalt[i,"region"] &
                                      bbsalt1$BCR == bbsalt[i,"BCR"] &
                                      bbsalt1$route == bbsalt[i,"route"] &
                                      bbsalt1$year == bbsalt[i,"year"]),"coni"],na.rm = F)
}
bbsalt$route = paste(bbsalt$BCR,bbsalt$route,sep = "_")


bbs = bbsalt[which(bbsalt$BCR %in% bcrs),]
maxbyrt = tapply(bbs$count,bbs$route,max,na.rm = T)
bbs = bbs[which(bbs$route %in% names(maxbyrt)[which(maxbyrt > 0)]),]
### remove bbs routes where species never observed



### compile data into a single dataframe with columns
#count, year, strat, route, program

bbs$program = "bbs"
nj$program = "nightjar"
bbs$programf = 1
nj$programf = 2

nprogs = 2







alldat = rbind(bbs,nj)


alldat = alldat[which(alldat$year > 2005),]






alldat$yearx = (alldat$year - min(alldat$year))+1

nyears = max(alldat$yearx)




alldat$strat = alldat$BCR#paste(alldat$prov,alldat$BCR,sep = "-")



alldat$stratx = as.integer(factor(alldat$strat)) 
ncounts = nrow(alldat)


nstrata = max(alldat$stratx)


nroutes = matrix(nrow = nprogs,
                 ncol = nstrata)

for(p in 1:nprogs){
  for(s in 1:nstrata){
ws = which(alldat$programf == p & alldat$stratx == s)
  alldat[ws,"routex"] = as.integer(factor(alldat[ws,"route"])) #these are unique within strata
  #alldat[ws,"obserx"] = as.integer(factor(alldat[ws,"obser"])) #these are unique within strata
nroutes[p,s] = max(alldat[ws,"routex"])
#nobservers[p,s] = max(alldat[ws,"obserx"])
}
}

nrtsbbs = NA
for(j in 1:ncol(nroutes)){
  
nrtsbbs[j] = max(nroutes[,j])

}

nyearspred = simparam[sj,"nyearspred"]


### trend to simulate
B.true = simparam[sj,"B.true"] #log-scale, multiplicative annual change that results in a 30% decline over 10 years


xstrata = c(5,8,9,10,12,14)#extra BCRs to expand the nightjar program

nstratax = nstrata+length(xstrata)
nroutesx = c(nrtsbbs,as.integer(bbsrtsbybcr[paste(xstrata)]))
fyear = unique(alldat[which(alldat$year == fyearnj),"yearx"])


jags.dat = list(ncounts = ncounts,
                nroutes = nroutes,
                nyears = nyears,
                nrtsbbs = nrtsbbs,
                #nobs = nobservers,
                nprogs = nprogs,
                nstrata = nstrata,
                count = alldat$count,
                year = alldat$yearx,
                strat = alldat$stratx,
                program = alldat$programf,
                route = alldat$routex,
                #obser = alldat$obserx,
                nyearspred = nyearspred,
                B = B.true,
                nstratax = nstratax,
                nroutesx = nroutesx,
                fyear = fyear)


# Bayesian model description -----------------------------------------------



newdir = paste0(getwd(),"/",sim)
dir.create(newdir)
#setwd(newdir)


mod.sh = "
model{
### priors
tau.year<- 1/pow(sd.year,2)# ~dgamma(0.01,0.01)## survey-wide year-effect precision 
sd.year ~ dunif(0.001,2) #<- 1/pow(tau.year,0.5)#

################
## #program specific intercepts
################


for(i in 1:nprogs){ 
prog[i] ~ dnorm(0,1)
###program specific observer variation
 tau.rte[i] ~ dgamma(0.001,0.001)
 sd.rte[i] <- 1/pow(tau.rte[i],0.5)

tau.noise[i] ~ dgamma(0.001,0.001)# <- 1/pow(sd.noise[i],2) ## 

sd.noise[i] <- 1/pow(tau.noise[i],0.5)#~ dunif(0.001,5)


################
## route-effects
################
for(s in 1:nstrata){ # the data are structured based on strata (bcrs) but the stratification is not used in the model
for(r in 1:nroutes[i,s]){
rte[i,s,r] ~ dnorm(prog[i],tau.rte[i])#
}

}#s

################
## extra route-effects
################
for(s in (nstrata+1):nstratax){
for(r in 1:nroutesx[s]){
rte[i,s,r] ~ dnorm(prog[i],tau.rte[i])#
}#rx

}#sx

}#i by prog

slope ~ dnorm(0,0.1) #long-term trend estimate

################
## year-effects
################
# yeareffect[1] ~ dnorm(0,0.01) # survey-wide 1st year
#  
# Optional first-difference year-effects approach used in preliminary models, generates similar results
# for(y in 2:nyears){
# yeareffect[y] ~ dnorm(yeareffect[y-1],tau.year)# survey-wide mean
# }# y

for(y in 1:nyears){
 yeareffect[y] ~ dnorm(0,tau.year)# random year effects - fluctuations around the trend-line
 }# y

################
## likelihood
################
for(i in 1:ncounts){
elambda1[i] <- (year[i]-fyear)*slope + yeareffect[year[i]] + rte[program[i],strat[i],route[i]] #
#elambda[i] ~ dt(elambda1[i],tau.noise, nu) #alternative heavy-tailed overdispersion
elambda[i] ~ dnorm(elambda1[i],tau.noise[program[i]]) #overdispersion
log(lambda[i]) <- elambda[i]
count[i] ~ dpois(lambda[i])
log_lik[i] <- logdensity.pois(count[i],lambda[i]) # log likelihood calculation used to calculate WAIC values
}

###################
### optional derived parameters
###################
# 
# for(y in 1:nyears){
# for(s in 1:nstrata){
# for(p in 1:nprogs){
# 
# for(r in 1:nroutes[p,s]){
# nr[s,y,p,r] <- exp(yeareffect[y] + rte[p,s,r] + 0.5*(1/tau.noise[p]))#prog[p] + strata[s] + 
# }
# 
# # for(r in 1:nrtsbbs[s]){
# # nryp[s,y,p,r] <- exp(yeareffect[y] + rte[1,s,r] + 0.5*(1/tau.noise[p]))#prog[p] + strata[s] + 
# # }
# 
# n[s,y,p] <- mean(nr[s,y,p,1:nroutes[p,s]])
# #nhyp[s,y,p] <- mean(nrhyp[s,y,p,1:nrtsbbs[s]])
# }
# 
# }
# }



##############################
### predicted timeseries of counts for next 11 years, with a 10% decline
##############################
# 
  
Btau <- 1/(pow((0.005),2)) #simulating route, year, strata, and program-specific trends centered on the rate required for a 30% decline

### sets true trend for each route and stratum, uses only nroutes for BBS because it has higher route numbers in all strata
for(s in 1:nstrata){

for(r in 1:nrtsbbs[s]){
betasM[s,r] ~ dnorm(B,Btau)
}

}

for(s in (1+nstrata):nstratax){

for(r in 1:nroutesx[s]){
betasM[s,r] ~ dnorm(B,Btau)
}

}


for(y in 1:nyearspred){
yeareffectpred[y] ~ dnorm(0,tau.year)

for(s in 1:nstrata){

for(p in 1:nprogs){
for(r in 1:nrtsbbs[s]){
betas[y,s,p,r] <- betasM[s,r]
}
for(r in 1:nroutes[p,s]){

elp1[s,y,p,r] <- ((betas[y,s,p,r]*y) + yeareffectpred[y] + rte[p,s,r] )#+ prog[p]+strata[s] + 
#elp2[s,y,p,r] ~ dt(elp1[s,y,p,r],tau.noise,nu) #heavy-tailed overdispersion
elp2[s,y,p,r] ~ dnorm(elp1[s,y,p,r],tau.noise[p]) #overdispersion
log(lambdap[s,y,p,r]) <- elp2[s,y,p,r]
countp[s,y,p,r] ~ dpois(lambdap[s,y,p,r])

}#r



}#p

}#s


for(s in (nstrata+1):nstratax){

for(p in 1:nprogs){
for(r in 1:nroutesx[s]){
betas[y,s,p,r] <- betasM[s,r]
elp1[s,y,p,r] <- ((betas[y,s,p,r]*y) + yeareffectpred[y] + rte[p,s,r] )#+ prog[p]+strata[s] + 
#elp2[s,y,p,r] ~ dt(elp1[s,y,p,r],tau.noise,nu) #heavy-tailed overdispersion
elp2[s,y,p,r] ~ dnorm(elp1[s,y,p,r],tau.noise[p]) #overdispersion
log(lambdap[s,y,p,r]) <- elp2[s,y,p,r]
countp[s,y,p,r] ~ dpois(lambdap[s,y,p,r])

}#r



}#p

}#s


}#y









}"


mod.name = "nightjar comined model w predictions new.txt"
cat(mod.sh,
    file = mod.name)

# jagsMod = jags.model( mod.name, 
#                       data= jags.dat ,  
#                       n.chains= 1 , 
#                       n.adapt= 500)


# Bayesian model fitting --------------------------------------------------


#MCMC settings
n.chain = 3
n.save = 10000
n.burn = 400000
n.thin = 500
n.iter = n.burn+(n.save*n.thin)

## parameters to monitor
params = c("slope",
           "sd.year",
           "sd.obs",
           "sd.year.st",
           "sd.noise",
           "sd.rte",
           "sd.obs",
           "prog",
           #"strata",
           "rte",
           "betasM")
t1 = Sys.time()

jagmod.comb = jags(data = jags.dat,
              n.chains = n.chain,
              n.iter = n.iter,
              n.thin = n.thin,
              n.burnin = n.burn,
              parameters.to.save = params,
              model.file = mod.name,
              parallel = T,
              modules = NULL)

t2 = Sys.time()

t2-t1

jagso.comb = as.data.frame(jagmod.comb$summary)


jagmod.comb.preds = update(jagmod.comb,
                          n.iter = 20000,
                          n.thin = 10,
                          parameters.to.save = "countp",
                          parallel = T,
                          codaOnly = "countp")











# 
pdf(file = paste0(newdir,"/jags summary combined w predictions simple.pdf"))
plot(jagmod.comb$samples)
dev.off()

pdf(file = paste0(newdir,"/jags summary combined w predictions post pred simple.pdf"))
plot(jagmod.comb.preds$samples)
dev.off()

write.csv(jagso.comb,paste0(newdir,"/jags summary combined model w predictions simple.csv"))



# save.image("combined model post jags w predictions decline.RData")
# 





# simulating trend estimates ----------------------------------------------



# load("combined model post jags w predictions decline.RData")
 
preds = jagmod.comb.preds$sims.list$countp
### simulate data from populations with 30% decline over 10 years
### preds[iter,strat,year,prog,route]
summary(preds[,,,1,],na.rm = T)
summary(preds[,,,2,],na.rm = T)

quantile(preds[,,,1,],na.rm = T,probs = c(0.001,0.01,0.025,0.05,seq(0.1,0.9,0.1),0.95,0.975,0.99,0.999))
quantile(preds[,,,2,],na.rm = T,probs = c(0.001,0.01,0.025,0.05,seq(0.1,0.9,0.1),0.95,0.975,0.99,0.999))

quantile(preds,na.rm = T,probs = c(0.001,0.01,0.025,0.05,seq(0.1,0.9,0.1),0.95,0.975,0.99,0.999))

overpreds = length(which(preds > (5*max(alldat$count))))
reasonablepreds = length(which(preds < (5*max(alldat$count))))
overpreds/reasonablepreds
#0.0004452661
## 0.04% of simulated observations are > 5 times larger than the largest observed count
# this line below removes these simulated counts treating them as missing data
# this is sub-optimal, but given that it affects simulated counts from both surveys,
# and it affects such a small percentage of the counts,
# it should have an important effect on the final results
preds[which(preds > (5*max(alldat$count)))] <- NA
### above removes counts that are more than 5X higher than the highest count ever observed on the nightjar survey
## the necessity of this highlights the inadequacies of the log-normal, extra-poisson overdispersion term
### with a sample as large as the one used in this simulation, just the variance in the coutns
### derived from the extra-poisson variation generates estimates of the mean count > 500 birds




simobs = melt(preds,
              varnames = c("iteration",
                           "strat",
                           "year",
                           "prog",
                           "rte"),
              na.rm = T,
              value.name = "simcount")

mnpred = tapply(simobs$simcount,simobs[,c("strat","year","prog")],mean,na.rm = T)
mdpred = tapply(simobs$simcount,simobs[,c("strat","year","prog")],quantile,probs = 0.5,na.rm = T)
lcipred = tapply(simobs$simcount,simobs[,c("strat","year","prog")],quantile,probs = 0.05,na.rm = T)
ucipred = tapply(simobs$simcount,simobs[,c("strat","year","prog")],quantile,probs = 0.95,na.rm = T)





out = data.frame(iteration = 1:dim(preds)[1])
pbar = progress_bar$new(total = dim(preds)[1])
for(i in 1:dim(preds)[1]){
  
  
  simdatt = preds[i,,,,]
  simdat1 = melt(simdatt,
                varnames = c("strat",
                             "year",
                             "prog",
                             "rte"),
                na.rm = T,
                value.name = "simcount")
  
  
  for(vers in c("","extra")){
    
    if(vers == ""){
      simdat = simdat1[which(simdat1$strat %in% 1:nstrata),]
    }else{
      simdat = simdat1
    }
    
  simdat$itr = factor(1:nrow(simdat))
  simdat$strat = factor(simdat$strat,ordered = F, levels = 1:max(simdat$strat))
  
resc = floor((nyearspred/2))
  
  #simdat$yr = (simdat$year-median(simdat$year))/10
  simdat$yr = (simdat$year-median(simdat$year))/resc
  
 
  
  
  simdatbbs = simdat[which(simdat$prog == 1),]
  simdatbbs$prog = factor(simdatbbs$prog)
  simdatbbs$rte = factor(simdatbbs$rte)
  
  
  
  simdatnj = simdat[which(simdat$prog == 2),]
  simdatnj$prog = factor(simdatnj$prog)
  simdatnj$rte = factor(simdatnj$rte)
  
  
  
  
  
  simdat$prog = factor(simdat$prog)
  simdat$rte = factor(simdat$rte)
  
  
  
  # tmbbs2 = try(glmer(formula = simcount~ yr*strat + (1|rte)+(1|itr),
  #                   data = simdatbbs,
  #                   family = poisson(link = "log")),silent = T)
  tmbbs = try(glmer(formula = simcount~ yr + (1|rte)+(1|itr),
                    data = simdatbbs,
                    family = poisson(link = "log")),silent = T)
  # tmnj2 = try(glmer(formula = simcount~ yr*strat + (1|rte)+(1|itr),
  #                  data = simdatnj,
  #                  family = poisson(link = "log")),silent = T)
  tmnj = try(glmer(formula = simcount~ yr + (1|rte)+(1|itr),
                   data = simdatnj,
                   family = poisson(link = "log")),silent = T)#,
                   #control = glmerControl(optimizer = c("bobyqa")),
                   #check.conv.grad = .makeCC(action = "ignore",tol = 1e-6)),silent = T)

  
  # tmnj = try(glmer.nb(formula = simcount~ yr + (1|rte)+(1|itr),
  #                  data = simdatnj),silent = T)
  #            # ,
             #       family = poisson(link = "log"),
             #       control = glmerControl(optimizer = c("bobyqa")),
             #       check.conv.grad = .makeCC(action = "ignore",tol = 1e-6)),silent = T)
             # 
  
  
  # 
#   tm2 = try(glmer(formula = simcount~ yr*strat + (1|rte)+(1|itr),
#                    data = simdat,
#                    family = poisson(link = "log")),silent = T)
  tm = try(glmer(formula = simcount~ yr + prog + (1|rte)+(1|itr),
                 data = simdat,
                 family = poisson(link = "log")),silent = T)

  
  # tmbbs = glmer(formula = simcount~yr+(yr|strat)+(1|rte)+(1|itr),
  #               data = simdatbbs,
  #               family = poisson(link = "log"))
  # tmnj = glmer(formula = simcount~yr+(yr|strat)+(1|rte)+(1|itr),
  #              data = simdatnj,
  #              family = poisson(link = "log"))
  # tmbbs = glmmadmb(formula = simcount~yr+(yr|strat)+(1|rte)+(1|itr),
  #               data = simdatbbs,
  #               family = "poisson", link = "log")
  # tmnj = glmmadmb(formula = simcount~yr+(yr|strat)+(1|rte)+(1|itr),
  #              data = simdatnj,
  #              family = "poisson", link = "log")
  # 
  # 
  tocom2 = ""
  tocom = ""
  iired2 = NA 
  
  
  if(!is.null(tmnj@optinfo$conv$lme4$messages) |
     !is.null(tm@optinfo$conv$lme4$messages) |
     !is.null(tmbbs@optinfo$conv$lme4$messages)){
    tocom = "try reduced model"
    tmnj = try(glmer(formula = simcount~ yr + (1|itr),
                     data = simdatnj,
                     family = poisson(link = "log")),silent = T)
    tm = try(glmer(formula = simcount~ yr + prog + (1|itr),
                   data = simdat,
                   family = poisson(link = "log")),silent = T)
    tmbbs = try(glmer(formula = simcount~ yr + (1|itr),
                      data = simdatbbs,
                      family = poisson(link = "log")),silent = T)
    
    if(!is.null(tmnj@optinfo$conv$lme4$messages) |
       !is.null(tm@optinfo$conv$lme4$messages) |
       !is.null(tmbbs@optinfo$conv$lme4$messages)){
    
      reduced2 = T
      iired2 = 1
      while(reduced2){
    
      dnj2 = dplyr::sample_n(simdatnj,round(nrow(simdatnj)*1),replace = T)
      dbbs2 = dplyr::sample_n(simdatbbs,round(nrow(simdatbbs)*1),replace = T)
      dcomb2 = dplyr::sample_n(simdat,round(nrow(simdat)*1),replace = T)
      tmnj = try(glmer(formula = simcount~ yr + (1|itr),
                       data = dnj2,
                       family = poisson(link = "log")),silent = T)
      tm = try(glmer(formula = simcount~ yr + prog + (1|itr),
                     data = dcomb2,
                     family = poisson(link = "log")),silent = T)
      tmbbs = try(glmer(formula = simcount~ yr + (1|itr),
                        data = dbbs2,
                        family = poisson(link = "log")),silent = T)
      tocom2 = "second reduced model"
      
      if(is.null(tmnj@optinfo$conv$lme4$messages) &
         is.null(tm@optinfo$conv$lme4$messages) &
         is.null(tmbbs@optinfo$conv$lme4$messages)){
      reduced2 = F
      }
      iired2 = iired2+1
      }
        
    }
    
  }
  
  
  
  
  if(class(tmnj) == 'glmerMod'){
    out[i,paste0(vers,"nj_y_coef")] = (1/resc)*summary(tmnj)$coefficients["yr","Estimate"]
    out[i,paste0(vers,"nj_y_coef_SE")] = (1/resc)*summary(tmnj)$coefficients["yr","Std. Error"]
    out[i,paste0(vers,"nj_y_coef_p")] = summary(tmnj)$coefficients["yr","Pr(>|z|)"]
    #regional 
    # if(class(tmnj2) == 'glmerMod'){
    # tt = summary(tmnj2)$coefficients
    # rrs = c("yr",paste0("yr:","strat",2:nstrata))
    # ttt = tt[rrs,c("Estimate","Std. Error")]
    # 
    # out[i,paste0("nj_y_coef_strat_",1)] = ttt[1,"Estimate"]
    # out[i,paste0("nj_y_coef_SE_strat_",1)] = ttt[1,"Std. Error"]
    # 
    # for(ss in 2:max(as.numeric(simdat$strat))){
    #   out[i,paste0("nj_y_coef_strat_",ss)] = ttt[ss,"Estimate"]+ttt[1,"Estimate"]
    #   out[i,paste0("nj_y_coef_SE_strat_",ss)] = ttt[ss,"Std. Error"]
    #   
    #   
    # }#ss
    # }
  }
  
  
  if(class(tmbbs) == 'glmerMod'){
    out[i,paste0(vers,"bbs_y_coef")] = (1/resc)*summary(tmbbs)$coefficients["yr","Estimate"]
    out[i,paste0(vers,"bbs_y_coef_SE")] = (1/resc)*summary(tmbbs)$coefficients["yr","Std. Error"]
    out[i,paste0(vers,"bbs_y_coef_p")] = summary(tmbbs)$coefficients["yr","Pr(>|z|)"]
    #regional 
#     if(class(tmbbs2) == 'glmerMod'){
#     tt = summary(tmbbs2)$coefficients
#     rrs = c("yr",paste0("yr:","strat",2:nstrata))
# ttt = tt[rrs,c("Estimate","Std. Error")]
#     
# out[i,paste0("bbs_y_coef_strat_",1)] = ttt[1,"Estimate"]
# out[i,paste0("bbs_y_coef_SE_strat_",1)] = ttt[1,"Std. Error"]
# 
#     for(ss in 2:max(as.numeric(simdat$strat))){
#       out[i,paste0("bbs_y_coef_strat_",ss)] = ttt[ss,"Estimate"]+ttt[1,"Estimate"]
#     out[i,paste0("bbs_y_coef_SE_strat_",ss)] = ttt[ss,"Std. Error"]
#     
#     
#   }#ss
# }
  }
  
  
  if(class(tm) == 'glmerMod'){
    out[i,paste0(vers,"comb_y_coef")] = (1/resc)*summary(tm)$coefficients["yr","Estimate"]
    out[i,paste0(vers,"comb_y_coef_SE")] = (1/resc)*summary(tm)$coefficients["yr","Std. Error"]
    out[i,paste0(vers,"comb_y_coef_p")] = summary(tm)$coefficients["yr","Pr(>|z|)"]
  }
  out[i,paste0(vers,"comment")] = tocom
  out[i,paste0(vers,"comment2")] = tocom2
  out[i,paste0(vers,"nrandomSamples")] = iired2
  
  }#vers
  
  rm(list = ls()[which(ls() %in% c("tmbbs2","tmnj2","tm2","tmbbs","tmnj","tm","tt","ttt"))])
  write.csv(out[c(1:i),],paste0(newdir,"/simulation summary results simple.csv"),row.names = F)
  
  pbar$tick()
  
}#i


ncnts.bbs = length(which(!is.na(preds[1,1:3,1,1,])))
ncnts.nj = length(which(!is.na(preds[1,1:3,1,2,])))

print(paste("assumes in the western study area",ncnts.bbs,"routes run annually from the BBS and ",ncnts.nj,"routes run annually from the nightjar survey"))

extrancnts.bbs = length(which(!is.na(preds[1,,1,1,])))
extrancnts.nj = length(which(!is.na(preds[1,,1,2,])))

print(paste("assumes nationally ",extrancnts.bbs,"routes run annually from the BBS and ",extrancnts.nj,"routes run annually from the nightjar survey"))









save.image(paste0(newdir,"/post sim summary.RData"))


}







########################## plotting
########################## plotting
########################## plotting
########################## plotting
########################## plotting
########################## plotting

# Plotting ----------------------------------------------------------------


library(ggplot2)
library(RColorBrewer)


simparam = data.frame(simparam = c("base","10pdecrease","stable","20yr"),
                      nyearspred = c(11,11,11,21),
                      B.true = c(log(((0.7)^(1/10))),
                                 log(((0.9)^(1/10))),
                                 0,
                                 log(((0.7)^(1/10)))),
                      null1 = c(0,
                                0,
                                log(((0.7)^(1/10))),
                                0),
                      null2 = c(0,
                                log(((0.7)^(1/10))),
                                log(((0.9)^(1/10))),
                                0))

#simparam = simparam[c(1:3),] #removes the 20-year simulation, which showed almost 100% power for all comparisons







for(sj in 1:nrow(simparam)){
  set.seed(2019)
  #setwd(paste0("C:/nightjar/nightjarfinal"))
  
  sim = simparam[sj,"simparam"]
  newdir = paste0(getwd(),"/",sim)
  
  #newdir = paste0(getwd(),"/",sim)
  #setwd(newdir)
  
  
  B.true = simparam[sj,"B.true"] #log-scale, multiplicative annual change that results in a 30% decline over 10 years
  


out = read.csv(paste0(newdir,"/simulation summary results simple.csv"),stringsAsFactors = F)

 # out[,"comb_y_coef"] = out[,"nj_y_coef"]
 # out[,"comb_y_coef_SE"] = out[,"nj_y_coef_SE"]
 # out[,"comb_y_coef_p"] = out[,"nj_y_coef_p"]
models = c("bbs","nj","comb")
names(models) = c("All-bird","Targeted","Combined")

pthresh = 0.95 #(alpha = 0.1, 2-sided comparison)

for(m in paste0(c("","extra"),rep(models,each = 2))){
    for(n in c("null1","null2")){
  out[,paste0(m,"_y_sig_LT_",n)] = F
  w.sig = which(out[,paste0(m,"_y_coef")]+ (qnorm(pthresh)*out[,paste0(m,"_y_coef_SE")]) < simparam[sj,n])
  out[w.sig,paste0(m,"_y_sig_LT_",n)] = T
  #out[which(is.na(out[,paste0(m,"_y_coef_p")])),paste0(m,"_y_sig_neg_trend")] = NA
  
  out[,paste0(m,"_y_sig_GT_",n)] = F
  w.sig = which(out[,paste0(m,"_y_coef")]- (qnorm(pthresh)*out[,paste0(m,"_y_coef_SE")]) > simparam[sj,n])
  out[w.sig,paste0(m,"_y_sig_GT_",n)] = T
  #out[which(is.na(out[,paste0(m,"_y_coef_p")])),paste0(m,"_y_sig_pos_trend")] = NA
  
  # out[,paste0(m,"_y_trend_incl_true")] = F
  # out[which((((out[,paste0(m,"_y_coef_SE")]*1.96) + out[,paste0(m,"_y_coef")]) > B.true) &
  #             ((out[,paste0(m,"_y_coef")] - (out[,paste0(m,"_y_coef_SE")]*1.96)) < B.true)),paste0(m,"_y_trend_incl_true")] = T
  # out[which(is.na(out[,paste0(m,"_y_coef_p")])),paste0(m,"_y_trend_incl_true")] = NA
  # 
  
  out[,paste0(m,"_y_bias")] = abs(out[,paste0(m,"_y_coef")] - B.true)
   
  # for(sj in 1:nstrata){
  #   
  #   out[,paste0(m,"_y_sig_neg_trend_strat_",sj)] = F
  #   out[which(out[,paste0(m,"_y_coef_strat_",sj)]+1.96*(out[,paste0(m,"_y_coef_SE_strat_",sj)]) < 0),paste0(m,"_y_sig_neg_trend_strat_",sj)] = T
  #   out[which(is.na(out[,paste0(m,"_y_coef_strat_",sj)])),paste0(m,"_y_sig_neg_trend_strat_",sj)] = NA
  #   
  #   out[,paste0(m,"_y_sig_pos_trend_strat_",sj)] = F
  #   out[which(out[,paste0(m,"_y_coef_strat_",sj)]-1.96*(out[,paste0(m,"_y_coef_SE_strat_",sj)]) > 0),paste0(m,"_y_sig_pos_trend_strat_",sj)] = T
  #   out[which(is.na(out[,paste0(m,"_y_coef_strat_",sj)])),paste0(m,"_y_sig_pos_trend_strat_",sj)] = NA
  #   
  #   # out[,paste0(m,"_y_trend_incl_true")] = F
  #   # out[which((((out[,paste0(m,"_y_coef_SE")]*1.96) + out[,paste0(m,"_y_coef")]) > B.true) &
  #   #             ((out[,paste0(m,"_y_coef")] - (out[,paste0(m,"_y_coef_SE")]*1.96)) < B.true)),paste0(m,"_y_trend_incl_true")] = T
  #   # out[which(is.na(out[,paste0(m,"_y_coef_p")])),paste0(m,"_y_trend_incl_true")] = NA
  #   # 
  #   
  #   out[,paste0(m,"_y_bias_strat_",sj)] = abs(out[,paste0(m,"_y_coef_strat_",sj)] - B.true)
  #   
  #   
  #   
  #   
  #   
  # }
    }#n
   
}#m





# 
# 
# 
# vars = c("coef","coef_SE","coef_p","bias","sig_neg_trend","sig_pos_trend","trend_incl_true")
# names(vars) = c("Trend (slope)","SE of Trend","p-value","Absolute bias in trend","significant negative trend","significant positive trend","trend includes truth")
# 



SEwidththreshhigh = (3.5/(2*qnorm(0.975)))
SEwidththreshmod = (6.7/(2*qnorm(0.975)))




# varcolrs = brewer.pal(3,"Dark2")
# names(varcolrs) = names(models)

varcolrs = rep(grey(0.2),3)
names(varcolrs) = names(models)


### consider ggplot versions
# pm = ggplot(data = out,aes(x = ))





######### plot the probability of trend < null1 and null2
for(gg in c("GT","LT")){
for(n in c("null1","null2")){
cls = paste0(paste0(rep(c("","extra"),times = 3),rep(c("bbs","nj","comb"),each = 2)),"_y_sig_",gg,"_",n)

tmp = data.frame(surv = rep(c("All Bird","Targeted","Combined"),each = 2),
                 Invest = rep(c("Current","National"),times = 3),
                 mean = NA)
for(j in 1:length(cls)){
  tmp[j,"mean"] = length(which(out[,cls[j]]))/nrow(out)
  
}

 tmp$surv = factor(tmp$surv,levels = c("All Bird","Targeted","Combined"),ordered = T)


pm = ggplot(data = tmp,aes(x = surv,y = mean,fill = Invest))+
  theme_classic()+
  geom_col(position = "dodge")+
  labs(x = "",y = paste("Probability of Trend Significantly",gg,signif(simparam[sj,n],2)))+
  scale_y_continuous(limits = c(0,1),expand = expand_scale(mult = c(0,0)))

pdf(paste0(newdir,paste("/Probability of Trend Significantly",gg,signif(simparam[sj,n],2)),".pdf"),
    height = 8,
    width = 4)
print(pm)
dev.off()

}
}




######### plot the SE of trend
    cls = paste0(paste0(rep(c("","extra"),times = 3),rep(c("bbs","nj","comb"),each = 2)),"_y_coef_SE")
    
    tmp = data.frame(surv = rep(c("All Bird","Targeted","Combined"),each = 2),
                     Invest = rep(c("Current","National"),times = 3),
                     mean = NA)
    for(j in 1:length(cls)){
      tmp[j,"mean"] = 100*(exp(mean(out[,cls[j]]))-1)
      tmp[j,"lci"] = 100*(exp(quantile(out[,cls[j]],0.025))-1)
      tmp[j,"uci"] = 100*(exp(quantile(out[,cls[j]],0.975))-1)
      tmp[j,"lqrt"] = 100*(exp(quantile(out[,cls[j]],0.75))-1)
      tmp[j,"uqrt"] = 100*(exp(quantile(out[,cls[j]],0.25))-1)
      tmp[j,"med"] = 100*(exp(quantile(out[,cls[j]],0.5))-1)
      
    }
    
    tmp$surv = factor(tmp$surv,levels = c("All Bird","Targeted","Combined"),ordered = T)
    
    
    pm = ggplot(data = tmp,aes(x = surv,y = med,colour = Invest))+
      theme_classic()+
      geom_hline(aes(yintercept = SEwidththreshhigh),colour = grey(0.5))+
      geom_hline(aes(yintercept = SEwidththreshmod),colour = grey(0.5))+
      geom_linerange(aes(ymin = lci,ymax = uci),size = 1,position = position_dodge(width = 0.5))+
      geom_pointrange(aes(ymin = lqrt,ymax = uqrt),size = 1.5,fatten = 2,position = position_dodge(width = 0.5))+
      #geom_point(position = position_dodge(width = 0.5))+
      labs(x = "",y = paste("SE of Trend"))+
      scale_y_continuous(limits = rev(c(max(tmp$uci)*1.1,0)),expand = expand_scale(mult = c(0,0)))
    
    
    
    
    
    pdf(paste0(newdir,paste("/SE of trend.pdf")),
        height = 8,
        width = 4)
    print(pm)
    dev.off()
    



    ######### plot the bias in trend
    cls = paste0(paste0(rep(c("","extra"),times = 3),rep(c("bbs","nj","comb"),each = 2)),"_y_bias")
    
    tmp = data.frame(surv = rep(c("All Bird","Targeted","Combined"),each = 2),
                     Invest = rep(c("Current","National"),times = 3),
                     mean = NA)
    for(j in 1:length(cls)){
      tmp[j,"mean"] = 100*(exp(mean(out[,cls[j]]))-1)
      tmp[j,"lci"] = 100*(exp(quantile(out[,cls[j]],0.025))-1)
      tmp[j,"uci"] = 100*(exp(quantile(out[,cls[j]],0.975))-1)
      tmp[j,"lqrt"] = 100*(exp(quantile(out[,cls[j]],0.75))-1)
      tmp[j,"uqrt"] = 100*(exp(quantile(out[,cls[j]],0.25))-1)
      tmp[j,"med"] = 100*(exp(quantile(out[,cls[j]],0.5))-1)
      
    }
    
    tmp$surv = factor(tmp$surv,levels = c("All Bird","Targeted","Combined"),ordered = T)
    
    
    pm = ggplot(data = tmp,aes(x = surv,y = med,colour = Invest))+
      theme_classic()+
      # geom_hline(aes(yintercept = SEwidththreshhigh),colour = grey(0.5))+
      # geom_hline(aes(yintercept = SEwidththreshmod),colour = grey(0.5))+
      geom_linerange(aes(ymin = lci,ymax = uci),size = 1,position = position_dodge(width = 0.5))+
      geom_pointrange(aes(ymin = lqrt,ymax = uqrt),size = 1.5,fatten = 2,position = position_dodge(width = 0.5))+
      #geom_point(position = position_dodge(width = 0.5))+
      labs(x = "",y = paste("Absolute Error of Trend"))+
      scale_y_continuous(limits = rev(c(max(tmp$uci)*1.1,0)),expand = expand_scale(mult = c(0,0)))
    
    
    
    
    
    pdf(paste0(newdir,paste("/Absolute Error of trend.pdf")),
        height = 8,
        width = 4)
    print(pm)
    dev.off()
    
    
    
    ######### alternative violin bias plot
    tmp = data.frame(surv = rep(c("All Bird","All Bird","Targeted","Targeted","Combined","Combined"),each = nrow(out)),
                     Invest = rep(rep(c("Current","National"),times = 3),each = nrow(out)),
                     AbsError = c(out[,cls[1]],out[,cls[2]],out[,cls[3]],out[,cls[4]],out[,cls[5]],out[,cls[6]]),
                     stringsAsFactors = F) 
    tmp$AbsError = 100*(exp(tmp$AbsError)-1)
    tmp$surv = factor(tmp$surv,levels = c("All Bird","Targeted","Combined"),ordered = T)
    

    pm = ggplot(data = tmp,aes(x = surv,y = AbsError,fill = Invest),borders = grey(0.5))+
      theme_classic()+
      # geom_hline(aes(yintercept = SEwidththreshhigh),colour = grey(0.5))+
      # geom_hline(aes(yintercept = SEwidththreshmod),colour = grey(0.5))+
      #geom_violin(aes(ymin = lci,ymax = uci),size = 1,position = position_dodge(width = 0.5))+
      geom_violin(scale = "count",
                  position = position_dodge(width = 0.8))+
      #geom_point(position = position_dodge(width = 0.5))+
      labs(x = "",y = paste("Absolute Error of Trend"))+
      scale_y_continuous(limits = rev(c(max(tmp$AbsError),0)),expand = expand_scale(mult = c(0,0)))
    
    
    
    
    
    pdf(paste0(newdir,paste("/Violin Absolute Error of trend.pdf")),
        height = 8,
        width = 4)
    print(pm)
    dev.off()
    




write.csv(out, paste0(newdir,"/Nightjar simulation results summary.csv"))




}























