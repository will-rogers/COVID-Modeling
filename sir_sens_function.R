#### SIR for sensitivity vs type of viral load ####

library(ggplot2)
library(dplyr)
library(data.table)
library(mc2d)
library(abind)

sir_sens <- function (sims, S, I, R, N, inf.days = 14, day.asym = 5, day.sym = 7, beta_vec = beta_vec,
                      delta.t=1, tests, ppn_sympt=ppn_sympt, 
                      contacts = 7, compliance = 1, care.seeking = 1,  
                      sensitivity = .8, times = 1, form = form) {
  form <- form
  inf.days <- inf.days
  
  dN_SI <- rbinom(n=sims,size=S,
                     prob=1-exp(-beta_vec*(apply(I[,day.asym:inf.days],1,sum))/(10000)*delta.t)) 
  dN_II <- I[,1:(inf.days-1)]
  dN_IR <- I[,inf.days]

  # update classes
  S. <- S - dN_SI 
  I. <- I
  I.[,1] <- I.[,1] + dN_SI
  I.[,2:inf.days] <- I.[,2:inf.days] + dN_II
  I.[,1:(inf.days-1)] <- I.[,1:(inf.days-1)] - dN_II
  I.[,inf.days] <- I.[,inf.days] - dN_IR
  R. <- R + dN_IR
  
  out <- cbind(S., I., R.) 

  mod.vec <- modify_sensitivity(form, inf.days)
  newSymptReportedTrue <- rbinom(sims, I.[,day.sym], ppn_sympt) # randomly draw symtomatic individuals
  newSymptReported <- floor(newSymptReportedTrue*care.seeking*(1-((1-sensitivity)^times))*mod.vec[day.sym])
  
  tested <- atests <- rmultinomial(sims,tests,out)
  
  for (i in 1:sims){
    for (j in 1:(inf.days+2)){
      if (j %in% 2:(inf.days+1)){
        tested[i,j] <- rbinom(1, atests[i,j], (1-((1-sensitivity)^times))*compliance*mod.vec[j-1])
      }
      if (j == 1 | j == (inf.days+2)) {
        tested[i,j] <- 0
      }
    }
  }
  
  sympt.isolate <- matrix(0,nr=sims,nc=inf.days+2) # storage for symptomatic cases to isolate
  sympt.isolate[,day.sym] <- newSymptReported
  
  atests.isolate <- tested
  
  contacts.contacted <- rmultinomial(sims,
                           rep(rpois(sims,contacts)*(newSymptReported + apply(atests.isolate, 1, sum)),sims),
                           matrix(rep(1,inf.days+2),nr=sims,nc=inf.days+2,byrow=T)*out)
  
  out <- pmax(out - sympt.isolate - atests.isolate - contacts.contacted, 0)
  
  out <- cbind(out, atests, tested, newSymptReportedTrue, 
               newSymptReported, contacts.contacted, atests.isolate, sympt.isolate
  )
  # store all states -- SIR states plus tested, reported, contacts
}
