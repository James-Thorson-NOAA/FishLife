get_r = function(vec=NULL, Linf=exp(vec['Loo']), K=exp(vec['K']),
  t0=-0.1, W_a=0.001, W_b=3.04, tm=exp(vec['tm']), dm=tm/4,
  maxage=ceiling(min(100,2*exp(vec['tmax']))), h=vec['h'], M=exp(vec['M']) ){

  age_vec = 0:maxage
  nages = length(age_vec)

  #Length-at-age
  length_a = Linf*(1-exp(-K*(age_vec-t0)))
  # Biomass-at-age
  biomass_a = W_a * length_a^W_b
  # Maturity-at-age
  maturity_a = 1/(1+exp(-(age_vec-tm)/dm))
  # mean Spawner weight at age
  SB_a = biomass_a * maturity_a

  # compute unfished Spawning biomass per recruit (SBR0)
  n0_a = rep(0,nages)
  for(aI in 1:nages){
    if(aI==1) n0_a[aI] = 1
    if(aI>=2) n0_a[aI] = n0_a[aI-1]*exp(-M)
    if(aI==nages) n0_a[aI] = n0_a[aI]/(1-exp(-M))
  }
  # Net reproductive rate
  NR = sum(n0_a * SB_a)

  # Unfished spawning biomass per recruit
  SBR0 = sum( n0_a * biomass_a * maturity_a )
  # Reproductive output Rs for bonyfish
  Rs = 4*h / ( SBR0*(1-h) )

  # Make Leslie matrix
  LeslieM_numbers = mat.or.vec(nages,nages)
  LeslieM_numbers[1,] = Rs * SB_a

  #fill rest of Matrix with Survival
  for(i  in 2:nages){
    LeslieM_numbers[i,(i-1)] = exp(-M)
  }

  # return intrinsic rate of population increase r and generation GT
  r_numbers = log( Re(eigen(LeslieM_numbers,only.values=TRUE)$values[1]) )
  generation_time = sum(age_vec * n0_a * SB_a) / NR
  Return = c("intrinsic_growth_rate"=r_numbers, "generation_time"=generation_time)

  # Generate Leslie matrix in biomass
  if( FALSE ){
    Bratio_aa = outer( biomass_a, biomass_a^-1 )
    LeslieM_biomass = LeslieM_numbers * Bratio_aa
    r_biomass = log( Re(eigen(LeslieM_biomass,only.values=TRUE)$values[1]) )
    Return = c(Return, "r_biomass"=r_biomass)
  }

  return(Return)
}
