#### SIR Step Function for Model objective 1 - identify role of LAMP specificity

# We assume a student model with on and off campus students with slightly different 
# contact structures, but no real difference in transmission. We don't do much here 
# with delays. Instead, we are focusing on how the testing sensitivity affects 
# some important parameters. Here, we assume that PCR could take as much at 2 days 
# to run, but LAMP is functionally instantaneous. 

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

sir_lamp <- function (sims, S.on, E.on, I1.on, I2.on,  R.on, N.on, newSympt1.on, newSympt2.on, beta_vec.on = beta_vec.on,
                      S.off,E.off,I1.off,I2.off,R.off,N.off,newSympt1.off, newSympt2.off, beta_vec.off = beta_vec.off, 
                      theta, gamma_I1I2, gamma_I2R, delta.t=1, tests, ppn_sympt=ppn_sympt, 
                      contacts.on = 7, contacts.off = 7, compliance = 1, care.seeking = 1,  
                      atest.wait.3,atest.wait.2,atest.wait.1,contact.wait.3,contact.wait.2,contact.wait.1,
                      test.scenario = c("2 Days","1 Day","No Delay"), sensitivity = .8, specificity = .99, times = 1) {

  dN_SE.on <- rbinom(n=sims,size=S.on,
                     prob=1-exp(-beta_vec.on*(I1.on+I2.on+I1.off+I2.off)/(N.on+N.off)*delta.t)) + rbinom(sims,1,.1) # add random introductions
  dN_EI1.on <- rbinom(n=sims,size=E.on,
                      prob=1-exp(-theta*delta.t))
  dN_I1I2.on <- rbinom(n=sims,size=I1.on,
                       prob=1-exp(-gamma_I1I2*delta.t))
  dN_I2R.on <- rbinom(n=sims,size=I2.on,
                      prob=1-exp(-gamma_I2R*delta.t))
  dN_SE.off <- rbinom(n=sims,size=S.off,
                      prob=1-exp(-beta_vec.off*(I1.on+I2.on+I1.off+I2.off)/(N.on+N.off)*delta.t)) + rbinom(sims,1,.1) # add random introductions
  dN_EI1.off <- rbinom(n=sims,size=E.off,
                       prob=1-exp(-theta*delta.t))
  dN_I1I2.off <- rbinom(n=sims,size=I1.off,
                        prob=1-exp(-gamma_I1I2*delta.t))
  dN_I2R.off <- rbinom(n=sims,size=I2.off,
                       prob=1-exp(-gamma_I2R*delta.t))
  # update classes
  S.on. <- S.on - dN_SE.on 
  E.on. <- E.on + dN_SE.on - dN_EI1.on 
  I1.on. <- I1.on + dN_EI1.on - dN_I1I2.on
  I2.on. <- I2.on + dN_I1I2.on - dN_I2R.on
  newSympt2.on <- newSympt1.on
  newSympt1.on <- dN_I1I2.on
  newSymptReportedTrue.on <- rbinom(sims,newSympt2.on,ppn_sympt) # randomly draw symtomatic individuals
  newSymptReported.on <- floor(newSymptReportedTrue.on*care.seeking*(1-((1-sensitivity)^times)))
  R.on. <- R.on + dN_I2R.on
  
  S.off. <- S.off - dN_SE.off 
  E.off. <- E.off + dN_SE.off - dN_EI1.off 
  I1.off. <- I1.off + dN_EI1.off - dN_I1I2.off
  I2.off. <- I2.off + dN_I1I2.off - dN_I2R.off
  newSympt2.off <- newSympt1.off
  newSympt1.off <- dN_I1I2.off
  newSymptReportedTrue.off <- rbinom(sims,newSympt2.off,ppn_sympt) # randomly draw symtomatic individuals
  newSymptReported.off <- floor(newSymptReportedTrue.off*care.seeking*(1-((1-sensitivity)^times)))
  R.off. <- R.off + dN_I2R.off
  out <- cbind( S.on.,  E.on.,  I1.on.,  I2.on., R.on., dN_I1I2.on, 
                S.off.,  E.off.,  I1.off.,  I2.off., R.off., dN_I1I2.off ) # assume that I1->I2 is when cases become detectable
  
  avail.tests <- tests #-(newSymptReported.on + newSymptReported.off) # we know testing is limited. Its is necessary to consider that
  # limited tests should be devoted to symptomatic cases first and THEN asymptomatic. This code takes the 
  # "available" # of tests and subtracts out the number of students needing a test to verify symptomology
  #avail.tests <- ifelse(avail.tests<0, 0, avail.tests) # if this number is less than zero, we are overdrawn for testing
  atests <- rmultinomial(sims,avail.tests,out[,c(1:5,7:11)])
  tested <- atests
  for (i in 1:sims){
    for (j in 1:10){
      if (j %in% c(3,4,8,9)){
        tested[i,j] <- rbinom(1, atests[i,j], 1-((1-sensitivity)^times)*compliance)
      }
      if (j %in% c(1,2,5,6,7,10)){
        tested[i,j] <- rbinom(1, atests[i,j], 1-((specificity)^times)*compliance)
      }
    }
  }
  
  sympt.isolate <- matrix(0,nr=sims,nc=10) # storage for symptomatic cases to isolate
  sympt.isolate[,c(4)] <- newSymptReported.on
  sympt.isolate[,c(9)] <- newSymptReported.off
  
  atests.isolate <- tested # holder for which tests will be positive that need to be isolated 
  # atests.isolate[,c(1,2,5,6,7,10)] <- 0 # set non-infected classes to 0
  atests.isolate <- floor(atests.isolate)
  
  atest.wait.3 <- sir_simple_step(atest.wait.2,sims,
                                  I1.on, I2.on, I1.off, I2.off, N.on, N.off,
                                  theta, gamma_I1I2, gamma_I2R,
                                  beta_vec.on, beta_vec.off)
  atest.wait.2 <- sir_simple_step(atest.wait.1,sims,
                                  I1.on, I2.on, I1.off, I2.off, N.on, N.off,
                                  theta, gamma_I1I2, gamma_I2R,
                                  beta_vec.on, beta_vec.off)
  atest.wait.1 <- atests.isolate
  
  if(test.scenario == "2 Days") {
    # randomly draw the contacts from the different classes
    tot.contacts.on <- rmultinomial(sims,
                                    rep(rpois(sims,contacts.on)*(newSymptReported.on + apply(atest.wait.3[,c(1:5)], 1, sum)),sims),
                                    matrix(c(1,1,1,1,1,1,1,1,1,1),nr=sims,nc=10,byrow=T)*out[,c(1:5, 7:11)])
    # contacts.on <- floor(tot.contacts.on*compliance)
    tot.contacts.off <- rmultinomial(sims,
                                     rep(rpois(sims,contacts.off)*(newSymptReported.off + apply(atest.wait.3[,c(6:10)], 1, sum)),sims),
                                     matrix(c(1,1,1,1,1,1,1,1,1,1),nr=sims,nc=10,byrow=T)*out[,c(1:5, 7:11)])
  }
  
  if(test.scenario == "1 Day") {
    # randomly draw the contacts from the different classes
    tot.contacts.on <- rmultinomial(sims,
                                    rep(rpois(sims,contacts.on)*(newSymptReported.on + apply(atest.wait.2[,c(1:5)], 1, sum)),sims),
                                    matrix(c(1,1,1,1,1,1,1,1,1,1),nr=sims,nc=10,byrow=T)*out[,c(1:5, 7:11)])
    # contacts.on <- floor(tot.contacts.on*compliance)
    tot.contacts.off <- rmultinomial(sims,
                                     rep(rpois(sims,contacts.off)*(newSymptReported.off + apply(atest.wait.2[,c(6:10)], 1, sum)),sims),
                                     matrix(c(1,1,1,1,1,1,1,1,1,1),nr=sims,nc=10,byrow=T)*out[,c(1:5, 7:11)])
  }
  
  if(test.scenario == "No Delay") {
    # randomly draw the contacts from the different classes
    tot.contacts.on <- rmultinomial(sims,
                                    rep(rpois(sims,contacts.on)*(newSymptReported.on + apply(atest.wait.1[,c(1:5)], 1, sum)),sims),
                                    matrix(c(1,1,1,1,1,1,1,1,1,1),nr=sims,nc=10,byrow=T)*out[,c(1:5, 7:11)])
    # contacts.on <- floor(tot.contacts.on*compliance)
    tot.contacts.off <- rmultinomial(sims,
                                     rep(rpois(sims,contacts.off)*(newSymptReported.off + apply(atest.wait.1[,c(6:10)], 1, sum)),sims),
                                     matrix(c(1,1,1,1,1,1,1,1,1,1),nr=sims,nc=10,byrow=T)*out[,c(1:5, 7:11)])
  }
  
  if(!test.scenario %in% c("2 Days","1 Day","No Delay")) {
    out <- 0
    print("Need correct delay interval")
  }
  
  tot.contacts <- tot.contacts.on + tot.contacts.off
  contacts <- tot.contacts
  for (i in 1:sims){
    for (j in 1:10){
      if (j %in% c(3,4,8,9)){
        contacts[i,j] <- rbinom(1, tot.contacts[i,j], 1-((1-sensitivity)^times)*compliance)
      }
      if (j %in% c(1,2,5,6,7,10)){
        contacts[i,j] <- rbinom(1, tot.contacts[i,j], 1-((specificity)^times)*compliance)
      }
    }
  }
  
  contact.wait.3 <- sir_simple_step(contact.wait.2,sims,
                                    I1.on, I2.on, I1.off, I2.off, N.on, N.off,
                                    theta, gamma_I1I2, gamma_I2R,
                                    beta_vec.on, beta_vec.off)
  contact.wait.2 <- sir_simple_step(contact.wait.1,sims,
                                    I1.on, I2.on, I1.off, I2.off, N.on, N.off,
                                    theta, gamma_I1I2, gamma_I2R,
                                    beta_vec.on, beta_vec.off)
  contact.wait.1 <- contacts

  if(test.scenario == "2 Days") {
    out[,c(1:5,7:11)] <- pmax(out[,c(1:5,7:11)] - sympt.isolate - (contact.wait.3) - (atest.wait.3),0)
  }
  
  if(test.scenario == "1 Day") {
    out[,c(1:5,7:11)] <- pmax(out[,c(1:5,7:11)] - sympt.isolate - (contact.wait.2) - (atest.wait.2),0)
  }
  
  if(test.scenario == "No Delay") {
    out[,c(1:5,7:11)] <- pmax(out[,c(1:5,7:11)] - sympt.isolate - (contact.wait.1) - (atest.wait.1),0)
  }
  
  if(!test.scenario %in% c("2 Days","1 Day","No Delay")) {
    out <- 0
    print("Need correct delay interval")
  }
  
  out <- cbind(out, atests, newSympt1.on, newSympt1.off, newSympt2.on, newSympt2.off, newSymptReported.on, newSymptReported.off,
               contacts, tot.contacts, avail.tests, atests.isolate,
               sympt.isolate, newSymptReportedTrue.on, newSymptReportedTrue.off, 
               atest.wait.3,atest.wait.2,atest.wait.1,
               contact.wait.3,contact.wait.2,contact.wait.1
  )
  # store all states -- SIR states plus tested, reported, contacts
}
