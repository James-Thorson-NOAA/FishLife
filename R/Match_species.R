
#' Match species
#'
#' Match \code{genus_species} to a given row of \code{ParentChild_gz}
#'
#' @param genus_species Character vector for Genus_species (separated by underscore character)
#' @inheritParams Calculate_ratio

#' @return integer of row numbers of \code{ParentChild_gz} matching \code{genus_species}

#' @export
Match_species = function( genus_species="Sebastes jordani", Database=FishLife::FishBase_and_RAM, ParentChild_gz=Database$ParentChild_gz ){
  # Match full taxonomy from fishbase
  genus_species = strsplit( tolower(genus_species), split=c(" ","_"))[[1]]
  Which = which( tolower(rfishbase::fishbase$Genus)==genus_species[1] & tolower(rfishbase::fishbase$Species)==genus_species[2] )
  if( length(Which)!=1 ) stop("Couldn't match input in fishbase")
  match_taxonomy = full_taxonomy = rfishbase::fishbase[Which,c("Class","Order","Family","Genus","Species")]

  # Match in database
  Count = 0
  Group = NA
  while( is.na(Group) ){
    Group = match( paste(tolower(match_taxonomy),collapse="_"), tolower(ParentChild_gz[,'ChildName']) )
    if( is.na(Group)){ match_taxonomy [ncol(match_taxonomy) - Count] <- "predictive" } 
    Count <- Count + 1
  }
  message( "Found match: ", paste(match_taxonomy,collapse=", ") )

  # Return match
  Return = list( "GroupNum"=Group, "match_taxonomy"=match_taxonomy, "full_taxonomy"=full_taxonomy)
  return( Return )
}
