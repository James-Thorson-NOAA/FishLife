
#' @export
Find_ancestors = function( child_num ){
  family_nums = child_num
  while(TRUE){
    if( is.na(ParentChild_gz[rev(family_nums)[1],'ParentRowNumber'])==TRUE ) break()
    family_nums = c(family_nums, ParentChild_gz[rev(family_nums)[1],'ParentRowNumber'])
  }
  return( family_nums )
}
