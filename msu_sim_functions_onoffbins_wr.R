
library(ggplot2)
library(dplyr)
library(data.table)
library(mc2d)
library(abind)


# SEIR Step Model from M. Ferrari Source Code -----------------------------
sir_simple_step <- function(waiting.group, sims,
                            I1.on, I2.on, I1.off, I2.off, N.on, N.off,
                            theta, gamma_I1I2, gamma_I2R,
                            beta_vec.on, beta_vec.off, delta.t = 1) {
  
  dN_SE.on <- rbinom(n=sims,size=waiting.group[,1],prob=1-exp(-beta_vec.on*(I1.on+I2.on+I1.off+I2.off)/(N.on+N.off)*delta.t)) # add random introductions
  dN_EI1.on <- rbinom(n=sims,size=waiting.group[,2],prob=1-exp(-theta*delta.t))
  dN_I1I2.on <- rbinom(n=sims,size=waiting.group[,3],prob=1-exp(-gamma_I1I2*delta.t))
  dN_I2R.on <- rbinom(n=sims,size=waiting.group[,4],prob=1-exp(-gamma_I2R*delta.t))
  dN_SE.off <- rbinom(n=sims,size=waiting.group[,6],prob=1-exp(-beta_vec.off*(I1.on+I2.on+I1.off+I2.off)/(N.on+N.off)*delta.t)) # add random introductions
  dN_EI1.off <- rbinom(n=sims,size=waiting.group[,7],prob=1-exp(-theta*delta.t))
  dN_I1I2.off <- rbinom(n=sims,size=waiting.group[,8],prob=1-exp(-gamma_I1I2*delta.t))
  dN_I2R.off <- rbinom(n=sims,size=waiting.group[,9],prob=1-exp(-gamma_I2R*delta.t))
  
  # update classes
  waiting.group[,1] <- waiting.group[,1] - dN_SE.on 
  waiting.group[,2] <- waiting.group[,2] + dN_SE.on - dN_EI1.on 
  waiting.group[,3] <- waiting.group[,3] + dN_EI1.on - dN_I1I2.on
  waiting.group[,4] <- waiting.group[,4] + dN_I1I2.on - dN_I2R.on
  waiting.group[,5] <- waiting.group[,5] + dN_I2R.on
  
  waiting.group[,6] <- waiting.group[,6] - dN_SE.off 
  waiting.group[,7] <- waiting.group[,7] + dN_SE.off - dN_EI1.off 
  waiting.group[,8] <- waiting.group[,8] + dN_EI1.off - dN_I1I2.off
  waiting.group[,9] <- waiting.group[,9] + dN_I1I2.off - dN_I2R.off
  waiting.group[,10] <- waiting.group[,10] + dN_I2R.off
  return(waiting.group)
}


sir_step <- function (sims, S.on, E.on, I1.on, I2.on,  R.on, N.on, newSympt1.on, newSympt2.on, beta_vec.on = beta_vec.on,
                      S.off,E.off,I1.off,I2.off,R.off,N.off,newSympt1.off, newSympt2.off, beta_vec.off = beta_vec.off, 
                      theta, gamma_I1I2, gamma_I2R, delta.t=1, tests, ppn_sympt=ppn_sympt, 
                      contacts.on = 7, contacts.off = 7, compliance = 1, care.seeking = 1, 
                      atest.wait.6, atest.wait.5,  atest.wait.4, atest.wait.3, atest.wait.2, atest.wait.1,
                      contact.wait.6, contact.wait.5, contact.wait.4, contact.wait.3, contact.wait.2, contact.wait.1,
                      test.scenario = c("5 Days", "4 Days","3 Days","2 Days","1 Day","No Delay")
                      ) {
  # adapted from Aaron King's code
  #sims - number of stochastic simulations
  #S vector of susceptibles, length=sims
  #E vector of exposed, length=sims
  #I1 vector of pre symptomatic infecteds, length=sims
  #I2 vector of possibly symptomatic infecteds, length=sims
  #R vector of recoverdeds, length=sims 
  #N vector of population size, length=sims
  #newSympt1  counter for new possibly symptomatic individuals, to allow 2 day health seeking delay
  #newSympt2  counter for new possibly symptomatic individuals, to allow 2 day health seeking delay
  #beta  transmission rate
  #theta rate from exposed to infected 
  #gamma_I1I2 rate from pre-symptomatic infecteds to possibly symptomatic
  #gamma_I2R rate from possibly symptomatic to recovered
  #delta.t time step length default is 1 day timestep
  #tests number of asymptomatic tests per day
  #ppn_sympt proportion of possibly symptomatic individuals who are symptomatic
  #contacts  average number of contaccts per individual
  #compliance is a proportion of confirmed asymptomatic tests or their contacts who either fail to 
  # show for testing or dont comply with isolation
  #care.seeking is the proportion of students who seek care once symptomatic
  # browser()
  dN_SE.on <- rbinom(n=sims,size=S.on,prob=1-exp(-beta_vec.on*(I1.on+I2.on+I1.off+I2.off)/(N.on+N.off)*delta.t)) + rbinom(sims,1,.1) # add random introductions
  dN_EI1.on <- rbinom(n=sims,size=E.on,prob=1-exp(-theta*delta.t))
  dN_I1I2.on <- rbinom(n=sims,size=I1.on,prob=1-exp(-gamma_I1I2*delta.t))
  dN_I2R.on <- rbinom(n=sims,size=I2.on,prob=1-exp(-gamma_I2R*delta.t))
  dN_SE.off <- rbinom(n=sims,size=S.off,prob=1-exp(-beta_vec.off*(I1.on+I2.on+I1.off+I2.off)/(N.on+N.off)*delta.t)) + rbinom(sims,1,.1) # add random introductions
  dN_EI1.off <- rbinom(n=sims,size=E.off,prob=1-exp(-theta*delta.t))
  dN_I1I2.off <- rbinom(n=sims,size=I1.off,prob=1-exp(-gamma_I1I2*delta.t))
  dN_I2R.off <- rbinom(n=sims,size=I2.off,prob=1-exp(-gamma_I2R*delta.t))
  # update classes
  S.on. <- S.on - dN_SE.on 
  E.on. <- E.on + dN_SE.on - dN_EI1.on 
  I1.on. <- I1.on + dN_EI1.on - dN_I1I2.on
  I2.on. <- I2.on + dN_I1I2.on - dN_I2R.on
  newSympt2.on <- newSympt1.on
  newSympt1.on <- dN_I1I2.on
  newSymptReportedTrue.on <- rbinom(sims,newSympt2.on,ppn_sympt) # randomly draw symtomatic individuals
  newSymptReported.on <- floor(newSymptReportedTrue.on*care.seeking)
  R.on. <- R.on + dN_I2R.on
  
  S.off. <- S.off - dN_SE.off 
  E.off. <- E.off + dN_SE.off - dN_EI1.off 
  I1.off. <- I1.off + dN_EI1.off - dN_I1I2.off
  I2.off. <- I2.off + dN_I1I2.off - dN_I2R.off
  newSympt2.off <- newSympt1.off
  newSympt1.off <- dN_I1I2.off
  newSymptReportedTrue.off <- rbinom(sims,newSympt2.off,ppn_sympt) # randomly draw symtomatic individuals
  newSymptReported.off <- floor(newSymptReportedTrue.off*care.seeking)
  R.off. <- R.off + dN_I2R.off
  out <- cbind( S.on.,  E.on.,  I1.on.,  I2.on., R.on., dN_I1I2.on, 
                S.off.,  E.off.,  I1.off.,  I2.off., R.off., dN_I1I2.off ) # assume that I1->I2 is when cases become detectable
  
  avail.tests <- tests-(newSymptReported.on + newSymptReported.off ) # we know testing is limited. Its is necessary to consider that
  # limited tests should be devoted to symptomatic cases first and THEN asymptomatic. This code takes the 
  # "available" # of tests and subtracts out the number of students needing a test to verify symptomology
  avail.tests <- ifelse(avail.tests<0, 0, avail.tests) # if this number is less than zero, we are overdrawn for testing
  atests <- rmultinomial(sims,avail.tests,out[,c(1:5,7:11)])
  
  sympt.isolate <- matrix(0,nr=sims,nc=10) # stoarge for symoptomatic cases to isolate
  sympt.isolate[,c(4)] <- newSymptReported.on
  sympt.isolate[,c(9)] <- newSymptReported.off
  
  atests.isolate <- atests # holder for which tests will be positive that need to be isolated 
  atests.isolate[,c(1,2,5,6,7,10)] <- 0 # set non-infected classes to 0
  atests.isolate <- floor(atests.isolate*compliance)
  
  atest.wait.6 <- sir_simple_step(atest.wait.5,sims,
                                  I1.on, I2.on, I1.off, I2.off, N.on, N.off,
                                  theta, gamma_I1I2, gamma_I2R,
                                  beta_vec.on, beta_vec.off)
  atest.wait.5 <- sir_simple_step(atest.wait.4,sims,
                                  I1.on, I2.on, I1.off, I2.off, N.on, N.off,
                                  theta, gamma_I1I2, gamma_I2R,
                                  beta_vec.on, beta_vec.off)
  atest.wait.4 <- sir_simple_step(atest.wait.3,sims,
                                  I1.on, I2.on, I1.off, I2.off, N.on, N.off,
                                  theta, gamma_I1I2, gamma_I2R,
                                  beta_vec.on, beta_vec.off)
  atest.wait.3 <- sir_simple_step(atest.wait.2,sims,
                                  I1.on, I2.on, I1.off, I2.off, N.on, N.off,
                                  theta, gamma_I1I2, gamma_I2R,
                                  beta_vec.on, beta_vec.off)
  atest.wait.2 <- sir_simple_step(atest.wait.1,sims,
                                  I1.on, I2.on, I1.off, I2.off, N.on, N.off,
                                  theta, gamma_I1I2, gamma_I2R,
                                  beta_vec.on, beta_vec.off)
  atest.wait.1 <- atests.isolate
  
  if(test.scenario == "5 Days") {
    # randomly draw the contacts from the different classes
    tot.contacts.on <- rmultinomial(sims,
                                    rep(rpois(sims,contacts.on)*(newSymptReported.on + apply(atest.wait.6[,c(3:4)], 1, sum)),sims),
                                    matrix(c(1,1,1,1,1,1,1,1,1,1),nr=sims,nc=10,byrow=T)*out[,c(1:5, 7:11)])
    # contacts.on <- floor(tot.contacts.on*compliance)
    tot.contacts.off <- rmultinomial(sims,
                                     rep(rpois(sims,contacts.off)*(newSymptReported.off + apply(atest.wait.6[,c(8:9)], 1, sum)),sims),
                                     matrix(c(1,1,1,1,1,1,1,1,1,1),nr=sims,nc=10,byrow=T)*out[,c(1:5, 7:11)])
  }
  
  if(test.scenario == "4 Days") {
    # randomly draw the contacts from the different classes
    tot.contacts.on <- rmultinomial(sims,
                                    rep(rpois(sims,contacts.on)*(newSymptReported.on + apply(atest.wait.5[,c(3:4)], 1, sum)),sims),
                                    matrix(c(1,1,1,1,1,1,1,1,1,1),nr=sims,nc=10,byrow=T)*out[,c(1:5, 7:11)])
    # contacts.on <- floor(tot.contacts.on*compliance)
    tot.contacts.off <- rmultinomial(sims,
                                     rep(rpois(sims,contacts.off)*(newSymptReported.off + apply(atest.wait.5[,c(8:9)], 1, sum)),sims),
                                     matrix(c(1,1,1,1,1,1,1,1,1,1),nr=sims,nc=10,byrow=T)*out[,c(1:5, 7:11)])
  }
  
  if(test.scenario == "3 Days") {
    # randomly draw the contacts from the different classes
    tot.contacts.on <- rmultinomial(sims,
                                    rep(rpois(sims,contacts.on)*(newSymptReported.on + apply(atest.wait.4[,c(3:4)], 1, sum)),sims),
                                    matrix(c(1,1,1,1,1,1,1,1,1,1),nr=sims,nc=10,byrow=T)*out[,c(1:5, 7:11)])
    # contacts.on <- floor(tot.contacts.on*compliance)
    tot.contacts.off <- rmultinomial(sims,
                                     rep(rpois(sims,contacts.off)*(newSymptReported.off + apply(atest.wait.4[,c(8:9)], 1, sum)),sims),
                                     matrix(c(1,1,1,1,1,1,1,1,1,1),nr=sims,nc=10,byrow=T)*out[,c(1:5, 7:11)])
  }
  
  if(test.scenario == "2 Days") {
    # randomly draw the contacts from the different classes
    tot.contacts.on <- rmultinomial(sims,
                                    rep(rpois(sims,contacts.on)*(newSymptReported.on + apply(atest.wait.3[,c(3:4)], 1, sum)),sims),
                                    matrix(c(1,1,1,1,1,1,1,1,1,1),nr=sims,nc=10,byrow=T)*out[,c(1:5, 7:11)])
    # contacts.on <- floor(tot.contacts.on*compliance)
    tot.contacts.off <- rmultinomial(sims,
                                     rep(rpois(sims,contacts.off)*(newSymptReported.off + apply(atest.wait.3[,c(8:9)], 1, sum)),sims),
                                     matrix(c(1,1,1,1,1,1,1,1,1,1),nr=sims,nc=10,byrow=T)*out[,c(1:5, 7:11)])
  }
  
  if(test.scenario == "1 Day") {
    # randomly draw the contacts from the different classes
    tot.contacts.on <- rmultinomial(sims,
                                    rep(rpois(sims,contacts.on)*(newSymptReported.on + apply(atest.wait.2[,c(3:4)], 1, sum)),sims),
                                    matrix(c(1,1,1,1,1,1,1,1,1,1),nr=sims,nc=10,byrow=T)*out[,c(1:5, 7:11)])
    # contacts.on <- floor(tot.contacts.on*compliance)
    tot.contacts.off <- rmultinomial(sims,
                                     rep(rpois(sims,contacts.off)*(newSymptReported.off + apply(atest.wait.2[,c(8:9)], 1, sum)),sims),
                                     matrix(c(1,1,1,1,1,1,1,1,1,1),nr=sims,nc=10,byrow=T)*out[,c(1:5, 7:11)])
  }
  
  if(test.scenario == "No Delay") {
    # randomly draw the contacts from the different classes
    tot.contacts.on <- rmultinomial(sims,
                                    rep(rpois(sims,contacts.on)*(newSymptReported.on + apply(atest.wait.1[,c(3:4)], 1, sum)),sims),
                                    matrix(c(1,1,1,1,1,1,1,1,1,1),nr=sims,nc=10,byrow=T)*out[,c(1:5, 7:11)])
    # contacts.on <- floor(tot.contacts.on*compliance)
    tot.contacts.off <- rmultinomial(sims,
                                     rep(rpois(sims,contacts.off)*(newSymptReported.off + apply(atest.wait.1[,c(8:9)], 1, sum)),sims),
                                     matrix(c(1,1,1,1,1,1,1,1,1,1),nr=sims,nc=10,byrow=T)*out[,c(1:5, 7:11)])
  }
  
  if(!test.scenario %in% c("5 Days", "4 Days","3 Days","2 Days","1 Day","No Delay")) {
    out <- 0
    print("Need correct delay interval")
  }
  tot.contacts <- tot.contacts.on + tot.contacts.off
  
  contacts <- tot.contacts
  contacts[,c(1,2,5,6,7,10)] <- 0
  contacts <- floor(contacts*compliance)
  
  contact.wait.6 <- sir_simple_step(contact.wait.5,sims,
                                    I1.on, I2.on, I1.off, I2.off, N.on, N.off,
                                    theta, gamma_I1I2, gamma_I2R,
                                    beta_vec.on, beta_vec.off)
  contact.wait.5 <- sir_simple_step(contact.wait.4,sims,
                                    I1.on, I2.on, I1.off, I2.off, N.on, N.off,
                                    theta, gamma_I1I2, gamma_I2R,
                                    beta_vec.on, beta_vec.off)
  contact.wait.4 <- sir_simple_step(contact.wait.3,sims,
                                    I1.on, I2.on, I1.off, I2.off, N.on, N.off,
                                    theta, gamma_I1I2, gamma_I2R,
                                    beta_vec.on, beta_vec.off)
  contact.wait.3 <- sir_simple_step(contact.wait.2,sims,
                                    I1.on, I2.on, I1.off, I2.off, N.on, N.off,
                                    theta, gamma_I1I2, gamma_I2R,
                                    beta_vec.on, beta_vec.off)
  contact.wait.2 <- sir_simple_step(contact.wait.1,sims,
                                    I1.on, I2.on, I1.off, I2.off, N.on, N.off,
                                    theta, gamma_I1I2, gamma_I2R,
                                    beta_vec.on, beta_vec.off)
  contact.wait.1 <- contacts
  # browser()
  if(test.scenario == "5 Days") {
    out[,c(1:5,7:11)] <- pmax(out[,c(1:5,7:11)] - sympt.isolate - (contact.wait.6) - (atest.wait.6),0)
  }
  
  if(test.scenario == "4 Days") {
    out[,c(1:5,7:11)] <- pmax(out[,c(1:5,7:11)] - sympt.isolate - (contact.wait.5) - (atest.wait.5),0)
  }

  if(test.scenario == "3 Days") {
    out[,c(1:5,7:11)] <- pmax(out[,c(1:5,7:11)] - sympt.isolate - (contact.wait.4) - (atest.wait.4),0)
  }

  if(test.scenario == "2 Days") {
    out[,c(1:5,7:11)] <- pmax(out[,c(1:5,7:11)] - sympt.isolate - (contact.wait.3) - (atest.wait.3),0)
  }

  if(test.scenario == "1 Day") {
    out[,c(1:5,7:11)] <- pmax(out[,c(1:5,7:11)] - sympt.isolate - (contact.wait.2) - (atest.wait.2),0)
  }
  
  if(test.scenario == "No Delay") {
    out[,c(1:5,7:11)] <- pmax(out[,c(1:5,7:11)] - sympt.isolate - (contact.wait.1) - (atest.wait.1),0)
  }
  
  if(!test.scenario %in% c("5 Days", "4 Days","3 Days","2 Days","1 Day","No Delay")) {
    out <- 0
    print("Need correct delay interval")
  }
  
  out <- cbind(out, atests, newSympt1.on, newSympt1.off, newSympt2.on, newSympt2.off, newSymptReported.on, newSymptReported.off,
               contacts, tot.contacts, avail.tests, atests.isolate,
               sympt.isolate, newSymptReportedTrue.on, newSymptReportedTrue.off, 
               atest.wait.6,atest.wait.5,atest.wait.4,atest.wait.3,atest.wait.2,atest.wait.1,
               contact.wait.6,contact.wait.5,contact.wait.4,contact.wait.3,contact.wait.2,contact.wait.1
               )
  # store all states -- SIR states plus tested, reported, contacts
}


############################################################################################################
############################################################################################################
############################################################################################################

msu_sim <- function(tst = 500, compliance = 0.75, introductions = 5, ppn_sympt = .8, 
                    care.seeking = 0.5, R0.on = 3, R0.off = 1.5, test.scenario = c("5 Days", "4 Days","3 Days","2 Days","1 Day","No Delay")){
  
  Tsim <- 100        # time to simulate over, we only care about start
  sims <- 100    # number of simulations
  ppn_sympt <- as.numeric(ppn_sympt)
  compliance <- as.numeric(compliance)
  care.seeking <- as.numeric(care.seeking) 
  introductions <- as.numeric(introductions)
  tst <- as.numeric(tst)
  test.scenario <- test.scenario
  R0.on <- as.numeric(R0.on)   
  R0.off <- as.numeric(R0.off)
  
  # storage vectors
  RE.on <- R0.on*.85  # RE assuming some fraction of population is already immune
  RE.off <- R0.off*.85
  beta.normal.on <- RE.on * (1/9)  # calculate Beta
  beta.normal.off <- RE.off * (1/9)  # calculate Beta
  # beta.outbreak <- RE * (1/9) * (1-distancing_reduction)   # calculate Beta
  beta_vec.on <- rep(beta.normal.on,sims)
  beta_vec.off <- rep(beta.normal.off,sims)
  theta <- 1/5  # 5 days from infection to infectious
  gamma_I1I2 <- 1/2 # 2 days asymptomatic infectious
  gamma_I2R <- 1/7 # 7 days infectious (this is probably too short)
  
  # storage 
  inf.on <- matrix(NA,Tsim,sims)  # storage for infectious class
  case.on <- matrix(NA,Tsim,sims) # storage for daily cases
  S.on <- matrix(round(4187*0.85),1,sims)  # start with 0.85 susceptible
  E.on <- matrix(floor(introductions*.25),1,sims)
  I1.on <- matrix(0,1,sims)                 # start with 5 asymptomatic infectious
  I2.on <- matrix(0,1,sims)
  R.on <- matrix(4187 - S.on,1,sims)
  sympt1.on <- matrix(0,1,sims)
  sympt2.on <- matrix(0,1,sims)
  symptrep.on <- matrix(0,1,sims)
  symptactual.on <- matrix(0,1,sims)
  new_contacts.on <- matrix(0,1,sims)
  isolation.on <- matrix(0,1,sims)
  quarantine.on <- matrix(0,1,sims)
  new_cases.on <- matrix(0,1,sims)
  true_test_positives.on <- matrix(0,1,sims)
  total_traces.on <- matrix(0,1,sims)
  comply_test_positives.on <- matrix(0,1,sims)
  N.on <- S.on+E.on+I1.on+I2.on+R.on
  
  inf.off <- matrix(NA,Tsim,sims)  # storage for infectious class
  case.off <- matrix(NA,Tsim,sims) # storage for daily cases
  S.off <- matrix(round(12563*0.85),1,sims)  # start with 0.85 susceptible
  E.off <- matrix(ceiling(introductions*.75),1,sims)
  I1.off <- matrix(0,1,sims)                 # start with 5 asymptomatic infectious
  I2.off <- matrix(0,1,sims)
  R.off <- matrix(12563 - S.off,1,sims)
  sympt1.off <- matrix(0,1,sims)
  sympt2.off <- matrix(0,1,sims)
  symptrep.off <- matrix(0,1,sims)
  symptactual.off <- matrix(0,1,sims)
  new_contacts.off <- matrix(0,1,sims)
  isolation.off <- matrix(0,1,sims)
  quarantine.off <- matrix(0,1,sims)
  new_cases.off <- matrix(0,1,sims)
  true_test_positives.off <- matrix(0,1,sims)
  total_traces.off <- matrix(0,1,sims)
  comply_test_positives.off <- matrix(0,1,sims)
  N.off <- S.off+E.off+I1.off+I2.off+R.off
  
  
  atest.wait.6 <- array(0,c(1,sims,10))
  atest.wait.5 <- array(0,c(1,sims,10))
  atest.wait.4 <- array(0,c(1,sims,10))
  atest.wait.3 <- array(0,c(1,sims,10))
  atest.wait.2 <- array(0,c(1,sims,10))
  atest.wait.1 <- array(0,c(1,sims,10))
  contact.wait.6 <- array(0,c(1,sims,10))
  contact.wait.5 <- array(0,c(1,sims,10))
  contact.wait.4 <- array(0,c(1,sims,10))
  contact.wait.3 <- array(0,c(1,sims,10))
  contact.wait.2 <- array(0,c(1,sims,10))
  contact.wait.1 <- array(0,c(1,sims,10))
  
  # distancing_reduction <- 0.5 # if trigger is crossed, NPIs are imposed and transmission is reduced by this fraction
  # qi_trigger <- numeric(sims)
  
  for(ts in 2:Tsim){
    # if (test.scenario == "FL") {
    #   if(ts <= 20) {
    #     tests <- tst*3
    #   }
    #   if(ts >20) {
    #     tests <- ((tst*Tsim)-(tst*3*20))/80
    #   }
    # }
    # if (test.scenario == "RL") {
    #   if(ts < 80) {
    #     tests <- ((tst*Tsim)-(tst*3*20))/80
    #   }
    #   if(ts >= 80) {
    #     tests <- tst*3
    #   }
    # }
    # if (test.scenario == "NP") {
    #   tests <- tst
    # }
    tests <- tst
    out <- sir_step(sims, 
                    S.on[ts-1,], E.on[ts-1,], I1.on[ts-1,], I2.on[ts-1,], R.on[ts-1,], 
                    N.on[ts-1,], sympt1.on[ts-1,], sympt2.on[ts-1,], beta_vec.on, 
                    S.off[ts-1,], E.off[ts-1,], I1.off[ts-1,], I2.off[ts-1,], R.off[ts-1,], 
                    N.off[ts-1,], sympt1.off[ts-1,], sympt2.off[ts-1,], beta_vec.off, 
                    theta, gamma_I1I2, gamma_I2R, delta.t=1, tests, contacts.on = 7, contacts.off = 7,
                    ppn_sympt = ppn_sympt, compliance = compliance, care.seeking = care.seeking,
                    atest.wait.6[ts-1,,],
                    atest.wait.5[ts-1,,],
                    atest.wait.4[ts-1,,],
                    atest.wait.3[ts-1,,],
                    atest.wait.2[ts-1,,],
                    atest.wait.1[ts-1,,],
                    contact.wait.6[ts-1,,],
                    contact.wait.5[ts-1,,],
                    contact.wait.4[ts-1,,],
                    contact.wait.3[ts-1,,],
                    contact.wait.2[ts-1,,],
                    contact.wait.1[ts-1,,],test.scenario) # call to SIR step function above
    S.on <- rbind(S.on,out[,1])  # update state
    E.on <- rbind(E.on,out[,2])  # update state
    I1.on <- rbind(I1.on,out[,3])  # update state
    I2.on <- rbind(I2.on,out[,4])  # update state
    R.on <- rbind(R.on,out[,5])  # update state
    N.on <- rbind(N.on,N.on[ts-1])  # update state
    sympt1.on <- rbind(sympt1.on,out[,23])  # update state
    sympt2.on <- rbind(sympt2.on,out[,25])  # update state
    symptrep.on <- rbind(symptrep.on,out[,27])  # update state
    symptactual.on <- rbind(symptactual.on,out[,70])
    new_contacts.on <- rbind(new_contacts.on,apply(out[,29:33],1,sum))  # update state
    true_test_positives.on <- rbind(true_test_positives.on, apply(out[,15:16],1,sum))  # update state
    comply_test_positives.on <- rbind(comply_test_positives.on, apply(out[,52:53],1,sum))  # update state
    new_cases.on <- rbind(new_cases.on,out[,6])  # update state
    total_traces.on <- rbind(total_traces.on, apply(out[,39:43],1,sum))# update state
    
    S.off <- rbind(S.off,out[,7])  # update state
    E.off <- rbind(E.off,out[,8])  # update state
    I1.off <- rbind(I1.off,out[,9])  # update state
    I2.off <- rbind(I2.off,out[,10])  # update state
    R.off <- rbind(R.off,out[,11])  # update state
    N.off <- rbind(N.off,N.off[ts-1])  # update state
    sympt1.off <- rbind(sympt1.off,out[,24])  # update state
    sympt2.off <- rbind(sympt2.off,out[,26])  # update state
    symptrep.off <- rbind(symptrep.off,out[,28])  # update state
    symptactual.off <- rbind(symptactual.off,out[,71])
    new_contacts.off <- rbind(new_contacts.off,apply(out[,34:38],1,sum))  # update state
    true_test_positives.off <- rbind(true_test_positives.off, apply(out[,20:21],1,sum))  # update state
    comply_test_positives.off <- rbind(comply_test_positives.off, apply(out[,57:58],1,sum))  # update state
    new_cases.off <- rbind(new_cases.off,out[,12])  # update state
    total_traces.off <- rbind(total_traces.off, apply(out[,44:48],1,sum))# update state
    
    atest.wait.6 <- abind(atest.wait.6, array(out[,72:81], c(1,sims,10)), along = 1)
    atest.wait.5 <- abind(atest.wait.5, array(out[,82:91], c(1,sims,10)), along = 1)
    atest.wait.4 <- abind(atest.wait.4, array(out[,92:101], c(1,sims,10)), along = 1)
    atest.wait.3 <- abind(atest.wait.3, array(out[,102:111], c(1,sims,10)), along = 1)
    atest.wait.2 <- abind(atest.wait.2, array(out[,112:121], c(1,sims,10)), along = 1)
    atest.wait.1 <- abind(atest.wait.1, array(out[,122:131], c(1,sims,10)), along = 1)
    contact.wait.6 <- abind(contact.wait.6, array(out[,132:141], c(1,sims,10)), along = 1)
    contact.wait.5 <- abind(contact.wait.5, array(out[,142:151], c(1,sims,10)), along = 1)
    contact.wait.4 <- abind(contact.wait.4, array(out[,152:161], c(1,sims,10)), along = 1)
    contact.wait.3 <- abind(contact.wait.3, array(out[,162:171], c(1,sims,10)), along = 1)
    contact.wait.2 <- abind(contact.wait.2, array(out[,172:181], c(1,sims,10)), along = 1)
    contact.wait.1 <- abind(contact.wait.1, array(out[,182:191], c(1,sims,10)), along = 1)
    ####################################################################################
    # Total in Isolation/Qurantine
    isolation.on <- rbind(isolation.on,apply(comply_test_positives.on[(max(1,ts-10)):ts,],2,sum) + apply(symptrep.on[(max(1,ts-10)):ts,],2,sum)) # isolate for 10 days
    quarantine.on <- rbind(quarantine.on, apply(new_contacts.on[(max(1,ts-14)):ts,],2,sum) ) # quarantine for 14 days
    
    isolation.off <- rbind(isolation.off,apply(comply_test_positives.off[(max(1,ts-10)):ts,],2,sum) + apply(symptrep.off[(max(1,ts-10)):ts,],2,sum)) # isolate for 10 days
    quarantine.off <- rbind(quarantine.off, apply(new_contacts.off[(max(1,ts-14)):ts,],2,sum) )
  }
  
  inf.on <- I1.on+I2.on   # total infectious on campus
  inf.off <- I1.off+I2.off   # total infectious on campus
  case.on <- new_cases.on # daily cases
  case.off <- new_cases.off
  
  return(list("active.inf.on" = inf.on, 
              "reporting.symptoms.on" = symptrep.on, 
              "all.symptomatics.on" = symptactual.on,
              "positive.asympt.on" = true_test_positives.on, 
              "new.cases.on" = case.on, 
              "isolation.complying.on" = isolation.on,
              "quarantine.complying.on" = quarantine.on,
              "total_traces.on" = total_traces.on,
              "active.inf.off" = inf.off, 
              "reporting.symptoms.off" = symptrep.off, 
              "all.symptomatics.off" = symptactual.off,
              "positive.asympt.off" = true_test_positives.off, 
              "new.cases.off" = case.off, 
              "isolation.complying.off" = isolation.off,
              "quarantine.complying.off" = quarantine.off,
              "total_traces.off" = total_traces.off,
              "tests"= matrix(tst,Tsim,sims),
              "compliance" = matrix(compliance,Tsim,sims),
              "introductions"= matrix(introductions,Tsim,sims),
              "ppn_sympt"= matrix(ppn_sympt,Tsim,sims),
              # 'thresh'= matrix(thresh,Tsim,sims),
              'care.seeking'= matrix(care.seeking,Tsim,sims),
              "R0.on" = matrix(R0.on,Tsim,sims),
              "R0.off" = matrix(R0.off,Tsim,sims),
              "test.scenario" = matrix(test.scenario,Tsim,sims)
              ))
}



msu_sims_par <- function(tests=seq(0,1500,500), compliance=c(1,.8,.65), introductions = 5, 
                         ppn_sympt = 0.8, care.seeking = 0.5, R0.on = 2.5, R0.off = 2.5, 
                         test.scenario = c("5 Days", "4 Days","3 Days","2 Days","1 Day","No Delay"), ncores=NULL){
  # browser()
  vars <- expand.grid('tests'=tests,
                      'compliance'=compliance,
                      'introductions'=introductions,
                      'ppn_sympt'=ppn_sympt,
                      'care.seeking'=care.seeking,
                      'R0.on'=R0.on,
                      'R0.off'=R0.off,
                      'test.scenario'=test.scenario)
  # vars <- vars %>% 
  #   filter(R0.on == R0.off)
  
  if (is.null(ncores)){
    output <- rbindlist(apply(vars,1,FUN=function(x) msu_sim(tst = x[1], compliance = x[2], introduction = x[3],
                                                             ppn_sympt = x[4], care.seeking = x[5],
                                                             R0.on = x[6], R0.off = x[7], test.scenario = x[8])))
  } else {
    cl <- parallel::makeCluster(ncores)
    parallel::clusterExport(cl,varlist=c('msu_sim','lengthen','sir_step'))
    parallel::clusterEvalQ(cl,{library(ggplot2)
                              library(dplyr)
                              library(data.table)
                              library(mc2d)})
    output <- rbindlist(parallel::parApply(cl,vars,1,
                                           FUN=function(x) msu_sim(tst = x[1], compliance = x[2], introduction = x[3],
                                                                   ppn_sympt = x[4], care.seeking = x[5],
                                                                   R0.on = x[6], R0.off = x[7], test.scenario = x[8])))
    parallel::stopCluster(cl)
    rm('cl')
  }
  output[,group:=rep(1:(.N/100),each=100)]
  output[,compliance:=paste(100*compliance,'%',sep='')]
  output[,tests:=tests]
  output[,ppn_sympt:=ppn_sympt]
  output[,R0.on:=R0.on]
  output[,R0.off:=R0.off]
  output[,day:=1:.N,by=group]
  return(output)
}


# Creating a running total for important metrics --------------------------

running_tots <- function(output) {
  return(output %>% 
    group_by(group, tests, compliance) %>% 
    mutate(cum_cases.on = cumsum(new.cases.on),
           cum_reporting.symptoms.on = cumsum(reporting.symptoms.on),
           cum_all.symptomatics.on = cumsum(all.symptomatics.on),
           cum_all.asymptomatics.on = cumsum(positive.asympt.on),
           cum_cases.off = cumsum(new.cases.off),
           cum_reporting.symptoms.off = cumsum(reporting.symptoms.off),
           cum_all.symptomatics.off = cumsum(all.symptomatics.off),
           cum_all.asymptomatics.off = cumsum(positive.asympt.off)
           )
    )
}

plotting.function <- function(output, var, log.10 = T) {
  if (log.10 == T) {
    plot(ggplot(output,aes_string(x = "day", y = var, color = "factor(tests)")) +
      geom_line(aes(group = group), alpha = 0.025) +
      geom_smooth() +
      facet_grid(compliance~.) +
      labs(
        x = "Day of School Year") +
      scale_y_log10(n.breaks = 10)+
      theme_classic())
  }
  if (log.10 == F) {
    plot(ggplot(output,aes_string(x = "day", y = var, color = "factor(tests)")) +
      geom_line(aes(group = group), alpha = 0.025) +
      geom_smooth() +
      facet_grid(compliance~.) +
      labs(
        x = "Day of School Year") +
      theme_classic())
  }
}


# summarize.output <- function(output, var, over.time = F) {
#     df_summ <- mtcars %>%
#       group_by_(.dots = group_var) %>%
#       summarise_(.dots = setNames(summ, summ_name))
#   if(over.time == F) {
#     summ.1 <- paste0('mean(', var, ')')
#     summ.2 <- paste0('sd(', var, ')')
#     summ_name.1 <- paste0('mean_', var)
#     summ_name.2 <- paste0('sd_', var)
#     out <- output %>% 
#       filter(day == max(day)) %>% 
#       group_by(.dots = c("tests", "compliance")) %>% 
#       summarize(.dots = setNames(summ.1, summ_name.1),
#                 .dots = setNames(summ.2, summ_name.2))
#   }
#   if(over.time == F) {
#     output %>% 
#     group_by(tests, compliance) %>% 
#       summarize(mean <- mean(as.name(var)),
#                 stdev <- sd(as.name(var)))
#   }
# }

plotting.summaries <- function (out) {
  out$lower <- out$mean - (1.96*out$sd)
  out$upper <- out$mean + (1.96*out$sd)
  out %>% 
    ggplot(aes(x=factor(tests), y = mean, fill = factor(tests))) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = lower, ymax = upper)) +
    facet_grid(.~compliance) + 
    scale_y_continuous(n.breaks = 10)+
    theme_classic()
}

