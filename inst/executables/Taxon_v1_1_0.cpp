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
matrix<Type> cov_matrix( vector<Type> L_val, int n_rows, int n_cols ){
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
  }
  // Diagonal matrix
  if( n_cols<=0 ){
    for(int r=0; r<n_rows; r++){
      Cov_rr(r,r) = L_val(Count)*L_val(Count);
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

  // Data
  DATA_MATRIX( Y_ij );
  DATA_IMATRIX( Missing_az );
  DATA_IMATRIX( PC_gz );
  DATA_IVECTOR( g_i );

  // Parameters
  PARAMETER_VECTOR( alpha_j );
  PARAMETER_VECTOR( L_z );
  PARAMETER_VECTOR( obsL_z );

  // Taxonomic random effects
  PARAMETER_MATRIX( beta_gj );

  // Missing data random effects
  PARAMETER_VECTOR( Y_a );

  // Derived data
  int n_j = Y_ij.row(0).size();
  int n_i = Y_ij.col(0).size();
  int n_g = PC_gz.col(0).size();

  // Objective funcction
  Type jnll = 0;
  vector<Type> jnll_comp( 6 );
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

  // Process covariance
  matrix<Type> Cov_jj( n_j, n_j );
  Cov_jj = cov_matrix(L_z, n_j, Options_vec(1));

  // Observation covariance
  matrix<Type> obsCov_jj( n_j, n_j );
  obsCov_jj = cov_matrix(obsL_z, n_j, Options_vec(0));

  // Probability of random effects
  vector<Type> Tmp_j( n_j );
  for( int g=0; g<n_g; g++ ){
    for( int j=0; j<n_j; j++ ){
      if( PC_gz(g,1)==0 ) Tmp_j(j) = beta_gj(g,j) - alpha_j(j);
      if( PC_gz(g,1)>=1 ) Tmp_j(j) = beta_gj(g,j) - beta_gj(PC_gz(g,0),j);
    }
    jnll_comp(PC_gz(g,1)) += MVNORM( Cov_jj )( Tmp_j );
  }

  // Probability of data
  matrix<Type> Yhat_ij( n_i, n_j );
  for( int i=0; i<n_i; i++){
    Yhat_ij.row( i ) = beta_gj.row( g_i(i) );
    jnll_comp(5) += MVNORM( obsCov_jj )( Ycomplete_ij.row(i) - Yhat_ij.row(i) );
  }

  // Return stuff
  REPORT( obsCov_jj );
  REPORT( Cov_jj );
  REPORT( beta_gj );
  REPORT( Yhat_ij );
  REPORT( Ycomplete_ij );
  REPORT( jnll );
  REPORT( jnll_comp );

  // Get covariances
  // ADREPORT( beta_gj );

  // Return jnll
  jnll = jnll_comp.sum();
  return jnll;
}
