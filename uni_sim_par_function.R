####### Function to iterate over a ton of different variables 


uni_sims_par <- function(tst = 500, 
                         compliance = 0.75, 
                         introductions = 5, 
                         ppn_sympt = .8, 
                         care.seeking = 0.5, 
                         R0.on = 3,  R0.off = 1.5, 
                         test.scenario = c("2 Days","1 Day","No Delay"),
                         sens.pcr = .99, spec.pcr = .99, 
                         sens.lamp = c(.8, .9, 1), spec.lamp = .99, 
                         lamp.diagnostic = F, 
                         size.intro.on = 1, prob.into.on =0.1,
                         size.intro.off = 1, prob.into.off =0.1,
                         immunity = 0.1, 
                         N0 = 16750, 
                         on.campus.prop = .25, 
                         contact.tracing.limit = 100,
                         pooling = 4, pooling.multi = 1,
                         days = 100, sims = 200,
                         engage.lamp = 1,
                         ncores=NULL){
  vars <- expand.grid("tst" = tst, 
                      "compliance" = compliance, 
                      "introductions" = introductions, 
                      "ppn_sympt" = ppn_sympt, 
                      "care.seeking" = care.seeking, 
                      "R0.on" = R0.on,  "R0.off" = R0.off, 
                      "test.scenario" = test.scenario,
                      "sens.pcr" = sens.pcr, "spec.pcr" = spec.pcr, 
                      "sens.lamp" = sens.lamp, "spec.lamp" = spec.lamp, 
                      "lamp.diagnostic" = lamp.diagnostic, 
                      "size.intro.on" = size.intro.on, "prob.into.on" =prob.into.on,
                      "size.intro.off" = size.intro.off, "prob.into.off" =prob.into.off,
                      "immunity" = immunity, 
                      "N0" = N0, 
                      "on.campus.prop" = on.campus.prop, 
                      "contact.tracing.limit" = contact.tracing.limit,
                      "pooling" = pooling, 
                      "pooling.multi" = pooling.multi,
                      "days" = days, "sims" = sims,
                      "engage.lamp" = engage.lamp)
  # vars <- vars %>% 
  #   filter(R0.on == R0.off)
  if (is.null(ncores)){
  output <- rbindlist(apply(vars,1,FUN=function(x) uni_sim(tst = x[1], 
                                                           compliance = x[2], 
                                                           introductions = x[3], 
                                                           ppn_sympt = x[4], 
                                                           care.seeking = x[5], 
                                                           R0.on = x[6],  R0.off = x[7], 
                                                           test.scenario = x[8],
                                                           sens.pcr = x[9], spec.pcr = x[10], 
                                                           sens.lamp = x[11], spec.lamp = x[12], 
                                                           lamp.diagnostic = x[13], 
                                                           size.intro.on = x[14], prob.into.on =x[15],
                                                           size.intro.off = x[16], prob.into.off =x[17],
                                                           immunity = x[18], 
                                                           N0 = x[19], 
                                                           on.campus.prop = x[20], 
                                                           contact.tracing.limit = x[21],
                                                           pooling = x[22],
                                                           pooling.multi = x[23],
                                                           days = x[24], sims = x[25],
                                                           engage.lamp = x[26])))
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
                                           FUN=function(x) uni_sim(tst = x[1], 
                                                                   compliance = x[2], 
                                                                   introductions = x[3], 
                                                                   ppn_sympt = x[4], 
                                                                   care.seeking = x[5], 
                                                                   R0.on = x[6],  R0.off = x[7], 
                                                                   test.scenario = x[8],
                                                                   sens.pcr = x[9], spec.pcr = x[10], 
                                                                   sens.lamp = x[11], spec.lamp = x[12], 
                                                                   lamp.diagnostic = x[13], 
                                                                   size.intro.on = x[14], prob.into.on =x[15],
                                                                   size.intro.off = x[16], prob.into.off =x[17],
                                                                   immunity = x[18], 
                                                                   N0 = x[19], 
                                                                   on.campus.prop = x[20], 
                                                                   contact.tracing.limit = x[21],
                                                                   pooling = x[22],
                                                                   pooling.multi = x[23],
                                                                   days = x[24], sims = x[25],
                                                                   engage.lamp = x[26])))
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



