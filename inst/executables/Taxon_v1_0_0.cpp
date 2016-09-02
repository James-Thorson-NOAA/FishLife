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

  // Data
  DATA_MATRIX( Z_ik );
  DATA_MATRIX( Y_ij );
  DATA_IMATRIX( Missing_az );

  // Parameters
  PARAMETER_VECTOR( alpha_j );
  //PARAMETER_ARRAY( L_lk );  // Each column is the elements of the (lower) cholesky decomposition of covariance for taxonomic level k
  PARAMETER_VECTOR( obsL_z );

  // Taxonomic random effects
  //PARAMETER_ARRAY( beta0_zj );
  //PARAMETER_ARRAY( beta1_zj );
  //PARAMETER_ARRAY( beta2_zj );
  //PARAMETER_ARRAY( beta3_zj );
  //PARAMETER_ARRAY( beta4_zj );

  // Missing data random effects
  PARAMETER_VECTOR( Y_a );

  // Derived data
  int n_k = Z_ik.row(0).size();
  int n_j = Y_ij.row(0).size();
  int n_i = Y_ij.col(0).size();

  // Objective funcction
  Type jnll = 0;
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

  // Covariances for phylogeny(Total: n_k)
  //array<Type> Cov_jjk( n_j, n_j, n_k );
  //for( int k=0; k<n_k; k++){
  //  Cov_jjk.col(k) = loadings_matrix(L_lk.col(k), n_k, n_k) * loadings_matrix(L_lk.col(k), n_k, n_k).transpose();
  //}

  // Observation covariance
  matrix<Type> obsCov_jj( n_j, n_j );
  obsCov_jj = cov_matrix(obsL_z, n_j, Options_vec(0));

  // Probability of random effects
  //for( int z=0; z<beta0_zj.col(0).size(); z++ ) jnll += MVNORM( Cov_jjk.col(0).matrix() )( beta0_zj.row(z) );
  //for( int z=0; z<beta1_zj.col(0).size(); z++ ) jnll += MVNORM( Cov_jjk.col(1).matrix() )( beta1_zj.row(z) );
  //for( int z=0; z<beta2_zj.col(0).size(); z++ ) jnll += MVNORM( Cov_jjk.col(2).matrix() )( beta2_zj.row(z) );
  //for( int z=0; z<beta3_zj.col(0).size(); z++ ) jnll += MVNORM( Cov_jjk.col(3).matrix() )( beta3_zj.row(z) );
  //for( int z=0; z<beta4_zj.col(0).size(); z++ ) jnll += MVNORM( Cov_jjk.col(4).matrix() )( beta4_zj.row(z) );

  // Probability of data
  matrix<Type> Yhat_ij( n_i, n_j );
  for( int i=0; i<n_i; i++){
    for( int j=0; j<n_j; j++ ){
      //Yhat_ij.row(i) = alpha_j + beta0_zj.row(Z_ik(i,0)) + beta1_zj.row(Z_ik(i,1)) + beta2_zj.row(Z_ik(i,2)) + beta3_zj.row(Z_ik(i,3)) + beta4_zj.row(Z_ik(i,4));
      Yhat_ij(i,j) = alpha_j(j);
    }
    jnll += MVNORM( obsCov_jj )( Ycomplete_ij.row(i) - Yhat_ij.row(i) );
  }

  // Return stuff
  REPORT( obsCov_jj );
  REPORT( Yhat_ij );
  REPORT( Ycomplete_ij );
  REPORT( jnll );

  return jnll;
}
