####### Function to iterate over a ton of different variables 


uni_sims_par <- function(tests=c(500), compliance=c(1), introductions = 25, 
                         ppn_sympt = 0.8, care.seeking = 1, R0.on = 2.5, R0.off = 2.5, 
                         test.scenario = c("2 Days","1 Day","No Delay"),
                         sens.pcr = .99, spec.pcr = .99, 
                         sens.lamp = c(.8, .9, 1), spec.lamp = .99, lamp.diagnostic = F, ncores=NULL){
  vars <- expand.grid('tests'=tests,
                      'compliance'=compliance,
                      'introductions'=introductions,
                      'ppn_sympt'=ppn_sympt,
                      'care.seeking'=care.seeking,
                      'R0.on'=R0.on,
                      'R0.off'=R0.off,
                      'test.scenario'=test.scenario,
                      "sens.pcr"=sens.pcr,
                      "spec.pcr"=spec.pcr,
                      "sens.lamp"=sens.lamp,
                      "spec.lamp"=spec.lamp,
                      "lamp.diagnostic"=lamp.diagnostic)
  # vars <- vars %>% 
  #   filter(R0.on == R0.off)
  if (is.null(ncores)){
  output <- rbindlist(apply(vars,1,FUN=function(x) uni_sim(tst = x[1], compliance = x[2], introduction = x[3],
                                                           ppn_sympt = x[4], care.seeking = x[5],
                                                           R0.on = x[6], R0.off = x[7], test.scenario = x[8],
                                                           sens.pcr = x[9], spec.pcr = x[10], sens.lamp = x[11], 
                                                           spec.lamp = x[12], lamp.diagnostic = x[13])))
  } else {
    cl <- parallel::makeCluster(ncores)
    parallel::clusterExport(cl,varlist=c('uni_sim','lengthen','sir_lamp', 'sir_simple_step'))
    parallel::clusterEvalQ(cl,{
      library(ggplot2)
      library(dplyr)
      library(data.table)
      library(mc2d)
      library(abind)
      source("uni_sim_function.R")})
    
    output <- rbindlist(parallel::parApply(cl,vars,1,
                                           FUN=function(x) uni_sim(tst = x[1], compliance = x[2], introduction = x[3],
                                                                   ppn_sympt = x[4], care.seeking = x[5],
                                                                   R0.on = x[6], R0.off = x[7], test.scenario = x[8],
                                                                   sens.pcr = x[9], spec.pcr = x[10], sens.lamp = x[11], 
                                                                   spec.lamp = x[12], lamp.diagnostic = x[13])))
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



