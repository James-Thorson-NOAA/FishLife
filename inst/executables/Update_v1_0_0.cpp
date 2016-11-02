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
  // Data
  DATA_MATRIX( Ynew_ij );
  DATA_IMATRIX( Missing_az );
  DATA_VECTOR( predMean_j );
  DATA_MATRIX( predCov_jj );
  DATA_MATRIX( obsCov_jj );

  // Updated effects
  PARAMETER_VECTOR( beta_j );

  // Missing data random effects
  PARAMETER_VECTOR( Y_a );

  // Derived data
  int n_j = Ynew_ij.row(0).size();
  int n_i = Ynew_ij.col(0).size();
  int n_a = Missing_az.col(0).size();

  // Objective funcction
  Type jnll = 0;
  vector<Type> jnll_comp( 2 );
  jnll_comp.setZero();
  using namespace density;

  // Complete data
  matrix<Type> Ycomplete_ij( n_i, n_j );
  Ycomplete_ij.setZero();
  for( int i=0; i<n_i; i++){
  for( int j=0; j<n_j; j++){
    if( !isNA(Ynew_ij(i,j)) ) Ycomplete_ij(i,j) = Ynew_ij(i,j);
  }}
  for( int a=0; a<n_a; a++){
    if( !isNA(Missing_az(a,0)) ){
      Ycomplete_ij( Missing_az(a,0), Missing_az(a,1) ) = Y_a(a);
    }
  }

  // Probability of predictive distribution
  jnll_comp(0) += MVNORM( predCov_jj )( beta_j - predMean_j );

  // Probability of data
  vector<Type> Tmp_j( n_j );
  for( int i=0; i<n_i; i++){
    for( int j=0; j<n_j; j++ ) Tmp_j(j) = beta_j(j) - Ycomplete_ij(i,j);
    jnll_comp(1) += MVNORM( obsCov_jj )( Tmp_j );
  }

  // Return stuff
  REPORT( beta_j );
  REPORT( Ycomplete_ij );

  // Return jnll
  jnll = jnll_comp.sum();
  return jnll;
}
