#### Setup sensivity vs viral load sims ####

sens_sim <- function(tst = 500, compliance = 0.75, introductions = 50, ppn_sympt = .8, 
                     care.seeking = 0.5, R0 = 3, sens.pcr = .99, spec.pcr = .99, 
                     sens.lamp = c(.8, .9, 1), spec.lamp = .99){
  
  Tsim <- 100        # time to simulate over, we only care about start
  sims <- 100    # number of simulations
  ppn_sympt <- as.numeric(ppn_sympt)
  compliance <- as.numeric(compliance)
  care.seeking <- as.numeric(care.seeking) 
  introductions <- as.numeric(introductions)
  tst <- as.numeric(tst)
  sensitivity <- as.numeric(sensitivity)
  times <- as.numeric(times)
  R0 <- as.numeric(R0)   
  form <- form
  inf.days = 14
  day.asym = 5
  day.sym = 7
  RE <- R0  # RE assuming some fraction of population is already immune
  beta.normal <- RE * (1/9)  # calculate Beta
  beta_vec <- rep(beta.normal,sims)
  
  # storage 
  inf <- matrix(NA,Tsim,sims)  # storage for infectious class
  case <- matrix(NA,Tsim,sims) # storage for daily cases
  S <- matrix(round(10000-introductions),1,sims)  # start with 0.85 susceptible
  I <- array(0,c(Tsim,sims,inf.days))
  I[1,,] <- rmultinomial(sims, introductions, rep(1, inf.days))
  R <- matrix(0,1,sims)
  symptrep <- matrix(0,1,sims)
  symptactual <- matrix(0,1,sims)
  new_contacts <- matrix(0,1,sims)
  active.inf <- matrix(0,1,sims)
  new_cases <- matrix(0,1,sims)
  true_test_positives <- matrix(0,1,sims)
  comply_test_positives <- matrix(0,1,sims)
  total_traces <- matrix(0,1,sims)
  comply_test_positives <- matrix(0,1,sims)
  new_cases <- matrix(0,1,sims)
  new_asym <- matrix(0,1,sims)
  new_sym <- matrix(0,1,sims)
  new_recov <- matrix(0,1,sims)
  isolation <- matrix(0,1,sims)
  quarantine <- matrix(0,1,sims)
  
  for(ts in 2:Tsim){
    tests <- tst
    out <- sir_sens(sims, S[ts-1,], I[ts-1,,], R[ts-1,], 
                    N[ts-1,], inf.days = inf.days, day.asym = day.asym, 
                    day.sym = day.sym, beta_vec = beta_vec,
                    delta.t=1, tests = tests, ppn_sympt=.2, 
                    contacts = 7, compliance = 1, care.seeking = 1,  
                    sensitivity = sensitivity,  times = times, form = form) # call to SIR step function above
    S <- rbind(S,out[,1])  # update state]
    I[ts,,] <- out[,2:(inf.days+1)]  # update state
    active.inf <- rbind(active.inf,apply(out[,2:(inf.days+1)],1,sum))
    R <- rbind(R,out[,inf.days+2])  # update state
    symptrep <- rbind(symptrep,out[,49])  # update state
    symptactual <- rbind(symptactual,out[,50])
    new_contacts <- rbind(new_contacts,apply(out[,51:66],1,sum))  # update state
    true_test_positives <- rbind(true_test_positives, apply(out[,17:32],1,sum))  # update state
    comply_test_positives <- rbind(comply_test_positives, apply(out[,33:48],1,sum))  # update state
    new_cases <- rbind(new_cases,out[,2])  # update state
    new_asym <- rbind(new_asym,out[,day.asym+1])
    new_sym <- rbind(new_asym,out[,day.sym+1])
    new_recov <- rbind(new_asym,out[,16])
    
    ####################################################################################
    # Total in Isolation/Qurantine
    isolation <- rbind(isolation,apply(comply_test_positives[(max(1,ts-10)):ts,],2,sum) + apply(symptrep[(max(1,ts-10)):ts,],2,sum)) # isolate for 10 days
    quarantine <- rbind(quarantine, apply(new_contacts[(max(1,ts-14)):ts,],2,sum) ) # quarantine for 14 days
  }

  return(list("S" = S,
              "I" = active.inf, 
              "R" = 10000-S-active.inf,
              "reporting.symptoms" = symptrep, 
              "all.symptomatics" = symptactual,
              "positive.asympt" = true_test_positives, 
              "new.cases" = new_cases, 
              "isolation.complying" = isolation,
              "quarantine.complying" = quarantine,
              "tests"= matrix(tst,Tsim,sims),
              "compliance" = matrix(compliance,Tsim,sims),
              "introductions"= matrix(introductions,Tsim,sims),
              "ppn_sympt"= matrix(ppn_sympt,Tsim,sims),
              'form'= matrix(form,Tsim,sims),
              'care.seeking'= matrix(care.seeking,Tsim,sims),
              "R0" = matrix(R0,Tsim,sims),
              "sensitivity" = matrix(sensitivity,Tsim,sims),
              "test.p.pat" = matrix(times,Tsim,sims),
              "inf.days" = inf.days,
              "day.asym" = day.asym,
              "day.sym" = day.sym
  ))
}
