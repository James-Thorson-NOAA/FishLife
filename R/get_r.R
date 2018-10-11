get_r = function(vec=NULL, Linf=exp(vec['Loo']), K=exp(vec['K']),
  t0=0, W_a=0.001, W_b=3.04, tm=exp(vec['tm']), dm=tm/4, minage=0,
  maxage=ceiling(min(100,2*exp(vec['tmax']))), h=vec['h'], M=exp(vec['M']) ){

  age_vec = minage:maxage
  nages = length(age_vec)

  #Length-at-age
  length_a = Linf*(1-exp(-K*(age_vec-t0)))
  # Weight-at-age
  weight_a = W_a * length_a^W_b
  # Maturity-at-age
  maturity_a = 1/(1+exp(-(age_vec-tm)/dm))
  # mean Spawner weight at age
  SB_a = weight_a * maturity_a

  # compute unfished Spawning biomass per recruit (SBR0)
  n0_a = rep(0,nages)
  for (t in 1:nages){
    if(t==1) n0_a[t] = 1
    if(t>1) n0_a[t] = n0_a[t-1]*exp(-M)
    if(t==nages) n0_a[t] = n0_a[t]/(1-exp(-M))
  }

  # Unfished spawning biomass per recruit
  SBR0 = sum( n0_a * weight_a * maturity_a )
  # Reproductive output Rs for bonyfish
  Rs = 4*h / ( SBR0*(1-h) )

  # Make Leslie matrix
  L.Mat=mat.or.vec(nages,nages)

  L.Mat[1,] = Rs * SB_a

  #fill rest of Matrix with Survival
  for(i  in 2:nages){
    L.Mat[i,(i-1)] = exp(-M)
  }
  # Net reproductive rate
  NR = sum(n0_a * SB_a)

  # return intrinsic rate of population increase r and generation GT
  intrinsic_growth_rate = log(as.numeric(eigen(L.Mat,only.values=TRUE)$values[1]))
  generation_time = sum(age_vec * n0_a * SB_a)/NR
  Return = c("intrinsic_growth_rate"=intrinsic_growth_rate, "generation_time"=generation_time)
  return(Return)
}
