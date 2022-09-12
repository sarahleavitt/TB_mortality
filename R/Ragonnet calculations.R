### calculate average mortality rate from Ragonnet et al, 2021
# this is from equation 2 in supplement
calcMortRate <- function(t,mort,cure,mu){
  
  denom <- cure+mort
  return(1-cure*exp(-mu*t)/denom-mort*exp(-cure*t-mort*t-mu*t)/denom)

}

# reported parameters in paper
mu <- 0.01 # background mortality rate without age
mort.sp  <- 0.389
mort.sn <- 0.025
cure.sp <- 0.231
cure.sn <- 0.130

mort.10.sn <- calcMortRate(10,mort.sn,cure.sn,mu)
mort.10.sp <- calcMortRate(10,mort.sp,cure.sp,mu)

mort.3.sn <- calcMortRate(3,mort.sn,cure.sn,mu)
mort.3.sp <- calcMortRate(3,mort.sp,cure.sp,mu)

### Calculate recovery at 3-4 years
# use equation 2 again
calcNatRec <- function(t,mort,cure,mu){
  
  a <- cure/(cure+mort)
  return(a*(exp(-mu*t)-exp(-t*(cure+mu+mort))))

}

cure.3.sp <- calcNatRec(3,mort.sp,cure.sp,mu)
cure.3.sn <- calcNatRec(3,mort.sn,cure.sn,mu)