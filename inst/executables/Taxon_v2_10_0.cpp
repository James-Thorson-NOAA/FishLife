#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Generate covariance matrix
// If positive, then factor decomposition
// If zero, then diagonal matrix
// If negative, then both factor decomposition and diagonal components
template<class Type>
matrix<Type> cov_matrix( vector<Type> L_val, vector<Type> logmult_col, Type min_var, int n_rows, int n_cols, bool invertTF ){
  // Define temporary objects
  matrix<Type> L_rc(n_rows, abs(n_cols));
  matrix<Type> Cov_rr(n_rows, n_rows);
  matrix<Type> Return_rr(n_rows, n_rows);
  Cov_rr.setZero();
  L_rc.setZero();
  // Loadings matrix with zero upper-diagonal
  int Count = 0;
  if( n_cols!=0 ){
    for(int r=0; r<n_rows; r++){
    for(int c=0; c<abs(n_cols); c++){
      if(r>=c){
        L_rc(r,c) = L_val(Count);
        Count++;
      }else{
        L_rc(r,c) = 0.0;
      }
    }}
    for(int c=0; c<abs(n_cols); c++){
      L_rc.col(c) = L_rc.col(c) * exp(logmult_col(c));
    }
  }
  // Diagonal matrix
  if( n_cols<=0 ){
    for(int r=0; r<n_rows; r++){
      Cov_rr(r,r) += L_val(Count)*L_val(Count);
      //Cov_rr(r,r) += exp( 2 * L_val(Count) );
      Count++;
    }
  }
  // Additive constant on diagonal
  for(int r=0; r<n_rows; r++){
    Cov_rr(r,r) += min_var;  // Necesary to prevent crashes during innner optimizer when using SR data
  }
  // Combine and return
  Cov_rr += L_rc * L_rc.transpose();
  if(invertTF==false) Return_rr = Cov_rr;
  if(invertTF==true) Return_rr = atomic::matinv( Cov_rr );
  return Return_rr;
}

// Main function
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Options
  DATA_IVECTOR( Options_vec );
  // Slot 0:  Number of observation factors
  // Slot 1:  Number of process error factors
  // Slot 2:  invertTF (whether to invert cov)
  // Slot 3:  Form for b_stock 0: bparam_stock = ln_b;  1: bparam_stock = log(phi_stock) =: log(SB_max / SB0)
  // Slot 4:  Turn off taxonomic hierarchy, such that Yhat_ij(i,j) = alpha_j(j) for all i and j
  DATA_VECTOR( Options );
  // Slot 0:  Additive constant for diagnonal of variance of obsCov_jj
  // Slot 1:  Additive constant for diagnonal of variance of Cov_jj
  // Slot 2:  SD_b_stock
  // Slot 3:  Penalty on low variance in ln_Rhat (to avoid B-H predictive function being flat over observed SSB range, which leads to singular inner hessian)
  DATA_IMATRIX( Cov_pz );

  // Data -- FishBase
  DATA_MATRIX( Y_ij );
  DATA_IMATRIX( Missing_az );
  DATA_IMATRIX( PC_gz );
  DATA_IVECTOR( g_i );

  // Data -- SR
  DATA_INTEGER( Nobs );
  DATA_INTEGER( Nstock );
  DATA_IVECTOR( Obs2Stock );
  DATA_IVECTOR( AR_Index );
  DATA_VECTOR( ln_R_obs );
  DATA_VECTOR( SSB_obs );
  DATA_VECTOR( SPRF0_stock );
  DATA_VECTOR( M_stock );
  DATA_VECTOR( SSBmax_stock );
  DATA_VECTOR( Rmax_stock );
  DATA_IVECTOR( j_SR );
  // j_SR(0):  Gives index j for ln_var
  // j_SR(1):  Gives index j for rho
  // j_SR(2):  Gives index j for ln_MASPS
  DATA_IVECTOR( i_stock );
  // i_stock(StockI):  Gives index i for stock Stock-Recruit stockI

  // Parameters -- FishBase
  PARAMETER_VECTOR( alpha_j );
  PARAMETER_VECTOR( L_z );
  PARAMETER_VECTOR( obsL_z );
  PARAMETER_VECTOR( L_logmult_col );
  PARAMETER_VECTOR( obsL_logmult_col );
  PARAMETER_VECTOR( cov_logmult_z ); // log-multiplier for process-error covariance for different taxonomic levels
  PARAMETER_MATRIX( beta_gj );
  PARAMETER_VECTOR( Y_a );

  // Parameters -- SR
  PARAMETER_VECTOR( bparam_stock );  // Nuissance parameter
  PARAMETER_VECTOR( gamma_p ); // Potential coefficients linking variables in beta_gj specified by Cov_design
  PARAMETER_VECTOR( theta_q ); // Potential coefficients linking variables in RAM database (turned off via mapping if absent)

  // Derived data
  int n_j = Y_ij.row(0).size();
  int n_i = Y_ij.col(0).size();
  int n_p = Cov_pz.col(0).size();
  int n_g = PC_gz.col(0).size();

  // Objective funcction
  Type jnll = 0;
  vector<Type> jnll_comp( 11 );
  jnll_comp.setZero();
  using namespace density;

  // Complete data
  matrix<Type> Ycomplete_ij( n_i, n_j );
  Ycomplete_ij.setZero();
  for( int i=0; i<n_i; i++){
  for( int j=0; j<n_j; j++){
    if( !isNA(Y_ij(i,j)) ) Ycomplete_ij(i,j) = Y_ij(i,j);
  }}
  for( int a=0; a<Missing_az.col(0).size(); a++){
    Ycomplete_ij( Missing_az(a,0), Missing_az(a,1) ) = Y_a(a);
  }

  /////////////////////////
  // FishBase part
  /////////////////////////

  // Process covariance
  matrix<Type> Cov_jj( n_j, n_j );
  Cov_jj = cov_matrix(L_z, L_logmult_col, Options(1), n_j, Options_vec(1), Options_vec(2));
  jnll_comp(8) = -1 * sum(dnorm( L_logmult_col, Type(0.0), Type(1.0), true ));

  // Observation covariance
  matrix<Type> obsCov_jj( n_j, n_j );
  obsCov_jj = cov_matrix(obsL_z, obsL_logmult_col, Options(0), n_j, Options_vec(0), Options_vec(2));
  jnll_comp(9) = -1 * sum(dnorm( obsL_logmult_col, Type(0.0), Type(1.0), true ));

  // Probability of random effects
  vector<Type> Parent_j( n_j );
  vector<Type> Prediction_j( n_j );
  vector<Type> Deviation_j( n_j );
  matrix<Type> tmpCov_jj( n_j, n_j );
  for( int g=0; g<n_g; g++ ){
    for( int j=0; j<n_j; j++ ){
      if( PC_gz(g,1)==0 ) Parent_j(j) = alpha_j(j);
      if( PC_gz(g,1)>=1 ) Parent_j(j) = beta_gj(PC_gz(g,0),j);
      Prediction_j(j) = Parent_j(j);
    }
    for( int p=0; p<n_p; p++ ){
      Prediction_j(Cov_pz(p,1)) += ( beta_gj(g,Cov_pz(p,0))-Parent_j(Cov_pz(p,0)) ) * gamma_p(p);
    }
    for( int j=0; j<n_j; j++ ){
      Deviation_j(j) = beta_gj(g,j) - Prediction_j(j);
    }
    tmpCov_jj = Cov_jj * exp(cov_logmult_z(PC_gz(g,1)));
    if( Options_vec(4)==false ){
      jnll_comp(PC_gz(g,1)) += MVNORM( tmpCov_jj )( Deviation_j );
    }
  }

  // Probability of data
  matrix<Type> Yhat_ij( n_i, n_j );
  for( int i=0; i<n_i; i++){
    if( Options_vec(4)==false ){
      Yhat_ij.row( i ) = beta_gj.row( g_i(i) );
    }else{
      Yhat_ij.row( i ) = alpha_j;
    }
    jnll_comp(5) += MVNORM( obsCov_jj )( Ycomplete_ij.row(i) - Yhat_ij.row(i) );
  }

  /////////////////////////
  // SR part
  /////////////////////////

  vector<Type> MLSPS_stock(Nstock);
  vector<Type> MASPS_stock(Nstock);
  vector<Type> h_stock(Nstock);
  vector<Type> logit_h_stock(Nstock);
  vector<Type> Fmsy_over_M_stock(Nstock);
  vector<Type> SPR_msy_stock(Nstock);
  vector<Type> ln_a_stock(Nstock);
  vector<Type> ro_stock(Nstock);
  vector<Type> SD_stock(Nstock);
  vector<Type> b_stock(Nstock);
  vector<Type> ln_M_stock(Nstock);
  ln_M_stock = log(M_stock);
  Type mean_ln_M_stock = ln_M_stock.sum() / ln_M_stock.size();

  vector<Type> mu_obs(Nobs);
  vector<Type> ln_R_obs_hat(Nobs);
  vector<Type> jnll_obs(Nobs);
  jnll_obs.setZero();

  // Extract parameter values from Ycomplete_ij
  for( int StockI=0; StockI<Nstock; StockI++ ){
    SD_stock(StockI) = exp( Ycomplete_ij(i_stock(StockI),j_SR(0)) / 2 );
    ro_stock(StockI) = Ycomplete_ij(i_stock(StockI),j_SR(1));
    MASPS_stock(StockI) = exp( Ycomplete_ij(i_stock(StockI),j_SR(2)) + theta_q(0) * (ln_M_stock(StockI)-mean_ln_M_stock) );
    MLSPS_stock(StockI) = MASPS_stock(StockI) / (1-exp(-M_stock(StockI)));
    // Ensure that MLSPS > 1, such that h > 0.2;  Changed from V2.5.0 to V2.6.0
    // MLSPS_stock(StockI) = MLSPS_stock(StockI) + 1;
    // Implies that Ycomplete_ij(StockI,j_SR(2)) is log-maximum lifetime spawners per spawner in excess of replacement, discounted to annual rate
    // Appears to result in failed inner-optimizer, perhaps because Ycomplete_ij(StockI,j_SR(2)) -> -Inf and gradient(Ycomplete_ij(StockI,j_SR(2))) -> 0 for some StockI
    h_stock(StockI) = MLSPS_stock(StockI) / ( 4 + MLSPS_stock(StockI) );
    logit_h_stock(StockI) = ( h_stock(StockI) - 0.2 ) / 0.8;
    Fmsy_over_M_stock(StockI) = pow( (4*h_stock(StockI)) / (1-h_stock(StockI)), 0.5 ) - 1;  // From Mangel et al. 2013 Eq. 13
    SPR_msy_stock(StockI) = pow( (1-h_stock(StockI)) / (4*h_stock(StockI)), 0.5 );          // From Mangel et al. 2013
    ln_a_stock(StockI) = log( MLSPS_stock(StockI) / SPRF0_stock(StockI) );   // From Myers et al. 1998, Eq. 5:  MLSPS=a*SPRF0, MASPS=MLSPS*(1-exp(-M)) ->  MASPS=a*SPRF0*(1-exp(-M)) -> a=MASPS/SPRF0/(1-exp(-M))
    if( j_SR.size()==3 ){
      if( Options_vec(3)==0 ){
        b_stock(StockI) = exp( bparam_stock(StockI) );
      }
      if( Options_vec(3)==1 ){
        // bparam_stock = R_max / R( SB -> Inf ) ->  beta = Rmax / alpha / bparam
        b_stock(StockI) = Rmax_stock(StockI) / exp(ln_a_stock(StockI)) / exp(bparam_stock(StockI));
      }
      if( Options_vec(3)==2 ){
        // bparam_stock = SB_max / SB(F=0)
        // crashes when exp(ln_a_stock(StockI)) * SPRF0_stock(StockI) < 1, i.e., MLSPS_stock < 1
        b_stock(StockI) = SSBmax_stock(StockI) / (exp(ln_a_stock(StockI)) * SPRF0_stock(StockI) - 1) / exp(bparam_stock(StockI));
      }
      if( Options_vec(3)==3 ){
        // bparam_stock _proportional-to_ SB_max
        b_stock(StockI) = SSBmax_stock(StockI) * exp(bparam_stock(StockI));
      }
      if( Options(2)>0 ) jnll_comp(6) -= dnorm( log(b_stock(StockI)), Type(0.0), Options(2), true );
    }
    if( j_SR.size()==4 ){
      if( Options_vec(3)==0 ){
        b_stock(StockI) = exp( Ycomplete_ij(i_stock(StockI),j_SR(3)) );
      }
      if( Options_vec(3)==1 ){
        // bparam_stock = R_max / R( SB -> Inf ) ->  beta = Rmax / alpha / bparam
        b_stock(StockI) = Rmax_stock(StockI) / exp(ln_a_stock(StockI)) / exp(Ycomplete_ij(i_stock(StockI),j_SR(3)));
      }
      if( Options_vec(3)==2 ){
        // bparam_stock = SB_max / SB(F=0)
        // crashes when exp(ln_a_stock(StockI)) * SPRF0_stock(StockI) < 1, i.e., MLSPS_stock < 1
        b_stock(StockI) = SSBmax_stock(StockI) / (exp(ln_a_stock(StockI)) * SPRF0_stock(StockI) - 1) / exp(Ycomplete_ij(i_stock(StockI),j_SR(3)));
      }
      if( Options_vec(3)==3 ){
        b_stock(StockI) = SSBmax_stock(StockI) * exp( Ycomplete_ij(i_stock(StockI),j_SR(3)) );
      }
    }
  }

  // Probability of data conditional on fixed and random effect values
  for( int ObsI=0; ObsI<Nobs; ObsI++ ){
    // Predict recruitment
    ln_R_obs_hat(ObsI) = ln_a_stock(Obs2Stock(ObsI)) + log( SSB_obs(ObsI) / ( 1 + SSB_obs(ObsI)/b_stock(Obs2Stock(ObsI)) ) );
    // Calculate divergence
    if( AR_Index(ObsI) == -1 ){
      mu_obs(ObsI) = ln_R_obs_hat(ObsI);
      jnll_obs(ObsI) = -1 * dnorm( ln_R_obs(ObsI), mu_obs(ObsI), SD_stock(Obs2Stock(ObsI)), true );
    }
    if( AR_Index(ObsI) == 0 ){
      mu_obs(ObsI) = ln_R_obs_hat(ObsI) + ro_stock(Obs2Stock(ObsI)) * ( ln_R_obs(ObsI-1) - ln_R_obs_hat(ObsI-1) );
      jnll_obs(ObsI) = -1 * dnorm( ln_R_obs(ObsI), mu_obs(ObsI), SD_stock(Obs2Stock(ObsI)), true );
    }
    if( AR_Index(ObsI) == 1 ){
      mu_obs(ObsI) = ln_R_obs_hat(ObsI);
      jnll_obs(ObsI) = -1 * dnorm( ln_R_obs(ObsI), mu_obs(ObsI), SD_stock(Obs2Stock(ObsI))/pow(1-pow(ro_stock(Obs2Stock(ObsI)),2),0.5), true );
    }
  }
  jnll_comp(7) = jnll_obs.sum();

  // Calculate CV in ln_R_obs_hat
  vector<Type> sum_ln_Rhat_stock(Nstock);
  vector<Type> mean_ln_Rhat_stock(Nstock);
  vector<Type> num_ln_Rhat_stock(Nstock);
  vector<Type> var_ln_Rhat_stock(Nstock);
  vector<Type> sd_ln_Rhat_stock(Nstock);
  sum_ln_Rhat_stock.setZero();
  num_ln_Rhat_stock.setZero();
  var_ln_Rhat_stock.setZero();
  for( int ObsI=0; ObsI<Nobs; ObsI++ ){
    num_ln_Rhat_stock(Obs2Stock(ObsI)) += 1;
    sum_ln_Rhat_stock(Obs2Stock(ObsI)) += ln_R_obs_hat(ObsI);
  }
  mean_ln_Rhat_stock = sum_ln_Rhat_stock / num_ln_Rhat_stock;
  for( int ObsI=0; ObsI<Nobs; ObsI++ ){
    var_ln_Rhat_stock(Obs2Stock(ObsI)) += pow( ln_R_obs_hat(ObsI) - mean_ln_Rhat_stock(Obs2Stock(ObsI)), 2 );
  }
  for( int StockI=0; StockI<Nstock; StockI++ ){
    sd_ln_Rhat_stock(StockI) = pow( var_ln_Rhat_stock(StockI) / num_ln_Rhat_stock(StockI), 0.5 );
  }

  // Penalize low variance in predictive recruitment
  jnll_comp(10) = -1 * abs(Options(3))*sum(log(sd_ln_Rhat_stock));
  REPORT( sum_ln_Rhat_stock );
  REPORT( mean_ln_Rhat_stock );
  REPORT( num_ln_Rhat_stock );
  REPORT( var_ln_Rhat_stock );
  REPORT( sd_ln_Rhat_stock );

  /////////////////////////
  // Reporting stuff
  /////////////////////////

  jnll = jnll_comp.sum();

  // Return stuff
  REPORT( obsCov_jj );
  REPORT( Cov_jj );
  REPORT( beta_gj );
  REPORT( Yhat_ij );
  REPORT( Ycomplete_ij );
  REPORT( jnll );
  REPORT( jnll_comp );
  REPORT( mu_obs );
  REPORT( ln_R_obs_hat );
  REPORT( jnll_obs );
  REPORT( ln_a_stock );
  REPORT( ln_R_obs_hat );
  REPORT( ro_stock );
  REPORT( SD_stock );
  REPORT( b_stock );

  // Get covariances
  // ADREPORT( beta_gj );
  if( Nstock>=2 ){
    REPORT( h_stock );
    REPORT( SPR_msy_stock );
    REPORT( Fmsy_over_M_stock );
    REPORT( logit_h_stock );
    // NOTE:  Using ADREPORT on any vector of length n_g causes memory problems!  message:  "Reached total allocation of 32673Mb: see help(memory.size)"
    //ADREPORT( SPR_msy );
    //ADREPORT( Fmsy_over_M );
    //ADREPORT( logit_h );
  }

  // Return jnll
  return jnll;
}
