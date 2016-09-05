
#' @export
Match_species = function( genus_species="Sebastes jordani", ParentChild_gz=Estimate_database$ParentChild_gz ){
  # Match full taxonomy from fishbase
  data( fishbase, package="rfishbase")
  genus_species = strsplit( tolower(genus_species), split=c(" ","_"))[[1]]
  Which = which( tolower(fishbase[,'Genus'])==genus_species[1] & tolower(fishbase[,'Species'])==genus_species[2] )
  if( length(Which)!=1 ) stop("Couldn't match input in fishbase")
  match_taxonomy = full_taxonomy = fishbase[Which,c("Class","Order","Family","Genus","Species")]

  # Match in database
  Count = 1
  Group = NA
  while( is.na(Group) ){
    Group = match( paste(tolower(match_taxonomy),collapse="_"), tolower(ParentChild_gz[,'ChildName']) )
    if( is.na(Group) ) rev(match_taxonomy)[Count] = "predictive"
  }
  message( "Found match: ", paste(match_taxonomy,collapse=", ") )

  # Return match
  Return = list( "GroupNum"=Group, "match_taxonomy"=match_taxonomy, "full_taxonomy"=full_taxonomy)
  return( Return )
}
