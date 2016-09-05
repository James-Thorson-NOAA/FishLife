
#' @export
Find_ancestors = function( child_num, ParentChild_gz=Estimate_database$ParentChild_gz ){

  # Search for all ancestors in the taxonomic tree
  family_nums = child_num
  while(TRUE){
    if( is.na(ParentChild_gz[rev(family_nums)[1],'ParentRowNumber'])==TRUE ) break()
    family_nums = c(family_nums, ParentChild_gz[rev(family_nums)[1],'ParentRowNumber'])
  }

  # Return vector of ancestors
  return( family_nums )
}
