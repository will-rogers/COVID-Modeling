#### University Simulation Functions


uni_sim <- function(tst = 500, compliance = 0.75, introductions = 5, ppn_sympt = .8, 
                    care.seeking = 0.5, R0.on = 3, R0.off = 1.5, test.scenario = c("2 Days","1 Day","No Delay"),
                    sens.pcr = .99, spec.pcr = .99, 
                    sens.lamp = c(.8, .9, 1), spec.lamp = .99, lamp.diagnostic = F){
  
  Tsim <- 100        # time to simulate over, we only care about start
  sims <- 200    # number of simulations
  lamp.diagnostic <- lamp.diagnostic
  ppn_sympt <- as.numeric(ppn_sympt)
  compliance <- as.numeric(compliance)
  care.seeking <- as.numeric(care.seeking) 
  introductions <- as.numeric(introductions)
  tst <- as.numeric(tst)
  test.scenario <- test.scenario
  sens.pcr <- as.numeric(sens.pcr)
  spec.pcr <- as.numeric(spec.pcr)
  sens.lamp <- as.numeric(sens.lamp)
  spec.lamp <- as.numeric(spec.lamp)
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
  S.on <- matrix(round(6287 * 0.85),1,sims)  # start with 0.85 susceptible
  E.on <- matrix(floor(introductions*.25),1,sims)
  I1.on <- matrix(0,1,sims)                 # start with 5 asymptomatic infectious
  I2.on <- matrix(0,1,sims)
  R.on <- matrix(6287 - S.on,1,sims)
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
  S.off <- matrix(round(10478 * 0.85),1,sims)  # start with 0.85 susceptible
  E.off <- matrix(ceiling(introductions*.75),1,sims)
  I1.off <- matrix(0,1,sims)                 # start with 5 asymptomatic infectious
  I2.off <- matrix(0,1,sims)
  R.off <- matrix(10478 - S.off,1,sims)
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
  symp.pcr <- matrix(0,1,sims)
  asymp.pcr <- matrix(0,1,sims)
  missed.pcr <- matrix(0,1,sims)
  N.off <- S.off+E.off+I1.off+I2.off+R.off
  
  atest.wait.3 <- array(0,c(1,sims,10))
  atest.wait.2 <- array(0,c(1,sims,10))
  atest.wait.1 <- array(0,c(1,sims,10))
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
    out <- sir_lamp(sims, 
                    S.on[ts-1,], E.on[ts-1,], I1.on[ts-1,], I2.on[ts-1,], R.on[ts-1,], 
                    N.on[ts-1,], sympt1.on[ts-1,], sympt2.on[ts-1,], beta_vec.on, 
                    S.off[ts-1,], E.off[ts-1,], I1.off[ts-1,], I2.off[ts-1,], R.off[ts-1,], 
                    N.off[ts-1,], sympt1.off[ts-1,], sympt2.off[ts-1,], beta_vec.off, 
                    theta, gamma_I1I2, gamma_I2R, delta.t=1, tests, contacts.on = 7, contacts.off = 7,
                    ppn_sympt = ppn_sympt, compliance = compliance, care.seeking = care.seeking,
                    atest.wait.3[ts-1,,],atest.wait.2[ts-1,,],atest.wait.1[ts-1,,],
                    contact.wait.3[ts-1,,], contact.wait.2[ts-1,,],contact.wait.1[ts-1,,],
                    test.scenario, sens.pcr = sens.pcr, spec.pcr = spec.pcr, 
                    sens.lamp = sens.lamp, spec.lamp = spec.lamp, lamp.diagnostic = lamp.diagnostic) # call to SIR step function above
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
    
    symp.pcr <- rbind(symp.pcr, out[,132])
    asymp.pcr <- rbind(asymp.pcr, out[,133])
    missed.pcr <- rbind(missed.pcr, out[,134])
    
    atest.wait.3 <- abind(atest.wait.3, array(out[,72:81], c(1,sims,10)), along = 1)
    atest.wait.2 <- abind(atest.wait.2, array(out[,82:91], c(1,sims,10)), along = 1)
    atest.wait.1 <- abind(atest.wait.1, array(out[,92:101], c(1,sims,10)), along = 1)
    contact.wait.3 <- abind(contact.wait.3, array(out[,102:111], c(1,sims,10)), along = 1)
    contact.wait.2 <- abind(contact.wait.2, array(out[,112:121], c(1,sims,10)), along = 1)
    contact.wait.1 <- abind(contact.wait.1, array(out[,122:131], c(1,sims,10)), along = 1)
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
              "symp.pcr" = symp.pcr,
              "asymp.pcr" = asymp.pcr,
              "missed.pcr" = missed.pcr,
              "tests"= matrix(tst,Tsim,sims),
              "compliance" = matrix(compliance,Tsim,sims),
              "introductions"= matrix(introductions,Tsim,sims),
              "ppn_sympt"= matrix(ppn_sympt,Tsim,sims),
              # 'thresh'= matrix(thresh,Tsim,sims),
              'care.seeking'= matrix(care.seeking,Tsim,sims),
              "R0.on" = matrix(R0.on,Tsim,sims),
              "R0.off" = matrix(R0.off,Tsim,sims),
              "test.scenario" = matrix(test.scenario,Tsim,sims),
              "sens.pcr" = matrix(sens.pcr,Tsim,sims),
              "spec.pcr" = matrix(spec.pcr,Tsim,sims),
              "sens.lamp" = matrix(sens.lamp,Tsim,sims),
              "spec.lamp" = matrix(spec.lamp,Tsim,sims),
              "lamp.diagnostic" = matrix(lamp.diagnostic,Tsim,sims)
              ))
}
