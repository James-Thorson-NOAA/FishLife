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
matrix<Type> cov_matrix( vector<Type> L_val, vector<Type> logmult_col, int n_rows, int n_cols ){
  // Define temporary objects
  matrix<Type> L_rc(n_rows, abs(n_cols));
  matrix<Type> Cov_rr(n_rows, n_rows);
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
      //Cov_rr(r,r) += L_val(Count)*L_val(Count);
      Cov_rr(r,r) += exp( 2 * L_val(Count) );
      Count++;
    }
  }
  // Combine and return
  Cov_rr += L_rc * L_rc.transpose();
  return Cov_rr;
}

// Main function
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Options
  DATA_IVECTOR( Options_vec );
  // Slot 0:  Number of observation factors
  // Slot 1:  Number of process error factors

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
  PARAMETER_VECTOR( ln_b_stock );  // Nuissance parameter

  // Derived data
  int n_j = Y_ij.row(0).size();
  int n_i = Y_ij.col(0).size();
  int n_g = PC_gz.col(0).size();

  // Objective funcction
  Type jnll = 0;
  vector<Type> jnll_comp( 10 );
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
  Cov_jj = cov_matrix(L_z, L_logmult_col, n_j, Options_vec(1));
  jnll_comp(8) = -1 * sum(dnorm( L_logmult_col, Type(0.0), Type(1.0), true ));

  // Observation covariance
  matrix<Type> obsCov_jj( n_j, n_j );
  obsCov_jj = cov_matrix(obsL_z, obsL_logmult_col, n_j, Options_vec(0));
  jnll_comp(9) = -1 * sum(dnorm( obsL_logmult_col, Type(0.0), Type(1.0), true ));

  // Probability of random effects
  vector<Type> Tmp_j( n_j );
  matrix<Type> tmpCov_jj( n_j, n_j );
  for( int g=0; g<n_g; g++ ){
    for( int j=0; j<n_j; j++ ){
      if( PC_gz(g,1)==0 ) Tmp_j(j) = beta_gj(g,j) - alpha_j(j);
      if( PC_gz(g,1)>=1 ) Tmp_j(j) = beta_gj(g,j) - beta_gj(PC_gz(g,0),j);
    }
    tmpCov_jj = Cov_jj * exp(cov_logmult_z(PC_gz(g,1)));
    jnll_comp(PC_gz(g,1)) += MVNORM( tmpCov_jj )( Tmp_j );
  }

  // Probability of data
  matrix<Type> Yhat_ij( n_i, n_j );
  for( int i=0; i<n_i; i++){
    Yhat_ij.row( i ) = beta_gj.row( g_i(i) );
    jnll_comp(5) += MVNORM( obsCov_jj )( Ycomplete_ij.row(i) - Yhat_ij.row(i) );
  }

  /////////////////////////
  // SR part
  /////////////////////////

  vector<Type> MLSPS(Nstock);
  vector<Type> ln_a_stock(Nstock);
  vector<Type> ro_stock(Nstock);
  vector<Type> SD_stock(Nstock);
  vector<Type> b_stock(Nstock);
  vector<Type> mu_obs(Nobs);
  vector<Type> ln_R_obs_hat(Nobs);
  vector<Type> jnll_obs(Nobs);
  jnll_obs.setZero();

  // Extract parameter values from Ycomplete_ij
  for( int StockI=0; StockI<Nstock; StockI++ ){
    SD_stock(StockI) = exp( Ycomplete_ij(i_stock(StockI),j_SR(0)) / 2 );
    ro_stock(StockI) = Ycomplete_ij(i_stock(StockI),j_SR(1));
    MLSPS(StockI) = exp( Ycomplete_ij(i_stock(StockI),j_SR(2)) ) / (1-exp(-M_stock(StockI)));
    ln_a_stock(StockI) = log( MLSPS(StockI) / SPRF0_stock(StockI) );   // From Myers et al. 1998, Eq. 5:  MLSPS=a*SPRF0, MASPS=MLSPS*(1-exp(-M)) ->  MASPS=a*SPRF0*(1-exp(-M)) -> a=MASPS/SPRF0/(1-exp(-M))
    b_stock(StockI) = exp( ln_b_stock(StockI) );
    jnll_comp(6) -= dnorm( ln_b_stock(StockI), Type(0.0), Type(10.0), true );
  }

  // Probability of data conditional on fixed and random effect values
  for( int ObsI=0; ObsI<Nobs; ObsI++ ){
    ln_R_obs_hat(ObsI) = ln_a_stock(Obs2Stock(ObsI)) + log( SSB_obs(ObsI) / ( 1 + SSB_obs(ObsI)/b_stock(Obs2Stock(ObsI)) ) );
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

  // Return jnll
  return jnll;
}
