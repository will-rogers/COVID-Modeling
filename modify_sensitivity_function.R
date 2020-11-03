#### Function of Viral Load probability modifiers on LAMP ####

modify_sensitivity <- function(form = c("uniform-max","uniform-mean","uniform-min",
                                        "normal", "early beta", "late beta"), inf.days = 14, range = c(0,1)){
  if(form == "uniform-max"){
    vect.mod <- rep(range[2], inf.days)
  }
  if(form == "uniform-mean"){
    vect.mod <- rep(mean(range), inf.days)
  }
  if(form == "uniform-min"){
    vect.mod <- rep(range[1], inf.days+1)
  }
  if(form == "normal"){
    vect.mod <- dnorm(1:14, inf.days/2, sd = inf.days/7)
    vect.mod <- vect.mod/max(vect.mod)
  }
  if(form == "early beta"){
    vect.mod <- dbeta(seq(0,1, length.out = inf.days), 2, 6)
    vect.mod <- vect.mod/max(vect.mod)
  }
  if(form == "late beta"){
    vect.mod <- dbeta(seq(0,1, length.out = inf.days), 3, 7)
    vect.mod <- vect.mod/max(vect.mod)
  }
  return(vect.mod)
}

a <- modify_sensitivity("uniform-max", 14)
b <- modify_sensitivity("normal", 14)
c <- modify_sensitivity("early beta", 14)
d <- modify_sensitivity("late beta", 14)

plot(1:14, a, type="l", col = "black", ylim = c(0,1))
lines(1:14,b, col = "blue")
lines(1:14,c, col = "red")
lines(1:14,d, col = "green")
abline(v=7)





