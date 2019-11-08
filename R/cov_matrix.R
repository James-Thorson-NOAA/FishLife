cov_matrix = function( L_val, n_rows, n_cols, output="Cov_rr" ){
  # Define temporary objects
  L_rc = matrix(nrow=n_rows, ncol=abs(n_cols))
  Cov_rr = matrix(nrow=n_rows, ncol=n_rows)
  Cov_rr[] = 0;
  L_rc[] = 0;
  # Loadings matrix with zero upper-diagonal
  Count = 1;
  if( n_cols!=0 ){
    for(r in 1:n_rows){
    for(c in 1:abs(n_cols)){
      if(r>=c){
        L_rc[r,c] = L_val[Count];
        Count = Count + 1;
      }else{
        L_rc[r,c] = 0.0;
      }
    }}
  }
  # Diagonal matrix
  if( n_cols<=0 ){
    for(r in 1:n_rows){
      Cov_rr[r,r] = Cov_rr[r,r] + L_val[Count]*L_val[Count];
      Count = Count + 1;
    }
  }
  # Combine and return
  Cov_rr = L_rc %*% t(L_rc);
  if(output=="Cov_rr") Return = Cov_rr
  if(output=="L_rc") Return = L_rc
  return( Return )
}
