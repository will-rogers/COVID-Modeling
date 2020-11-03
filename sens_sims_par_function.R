####### Function to iterate over a ton of different variables 

sens_sims_par <- function(tests=c(1,500,1000), compliance=c(1), introductions = 25, 
                         ppn_sympt = 0.8, care.seeking = 1, R0, sensitivity = c(0.8, 0.9, 1), 
                         times = c(1), form = c("uniform-max","uniform-mean","uniform-min",
                                                "normal", "early beta", "late beta"), ncores=NULL){
  vars <- expand.grid('tests'=tests,
                      'compliance'=compliance,
                      'introductions'=introductions,
                      'ppn_sympt'=ppn_sympt,
                      'care.seeking'=care.seeking,
                      'R0'=R0,
                      "sensitivity"=sensitivity,
                      "times"=times,
                      'form'=form)
  
  if (is.null(ncores)){
    output <- rbindlist(apply(vars,1,FUN=function(x) sens_sim(tst = x[1], compliance = x[2], introduction = x[3],
                                                              ppn_sympt = x[4], care.seeking = x[5],
                                                              R0 = x[6], sensitivity = x[7], times = x[8], form = x[9])))
  } else {
    cl <- parallel::makeCluster(ncores)
    parallel::clusterExport(cl,varlist=c('sens_sim', 'modify_sensitivity', 'sir_sens'))
    parallel::clusterEvalQ(cl,{
      library(ggplot2)
      library(dplyr)
      library(data.table)
      library(mc2d)
      library(abind)
      source("sens_sim_function.R")
      source("modify_sensitivity_function.R")})
    
    output <- rbindlist(parallel::parApply(cl,vars,1,
                                           FUN=function(x) sens_sim(tst = x[1], compliance = x[2], introduction = x[3],
                                                                    ppn_sympt = x[4], care.seeking = x[5],
                                                                    R0 = x[6], sensitivity = x[7], times = x[8], form = x[9])))
    parallel::stopCluster(cl)
    rm('cl')
  }
  output[,group:=rep(1:(.N/100),each=100)]
  output[,compliance:=paste(100*compliance,'%',sep='')]
  output[,tests:=tests]
  output[,ppn_sympt:=ppn_sympt]
  output[,R0:=R0]
  output[,day:=1:.N,by=group]
  return(output)
}
