#' Search species
#'
#' Match taxonomic inputs to a given row of \code{ParentChild_gz} or its closest ancestor
#'
#' This function attempts to do a smart match to elements of \code{ParentChild_gz}.  It sweeps from Order to Species
#' and ignores any taxonomic input listed as \code{"predictive"} until it finds something else.  It then appends
#' \code{"predictive"} to any lower taxonomic level that is missing, and checks whether this specification yields a single,
#' unique taxon.  If it does, it then returns the row number and potentially any ancestors (higher taxonomic levels)
#'
#' @param Class Character input for taxonomic class
#' @param Order Character input for taxonomic class
#' @param Family Character input for taxonomic class
#' @param Genus Character input for taxonomic class
#' @param Species Character input for taxonomic class
#' @param add_ancestors Boolean whether to add ancestors for matching species or not
#' @param ParentChild_gz vector providing index of parent-taxon for every child-taxa
#'
#' @return integer of row numbers of \code{ParentChild_gz} matching \code{genus_species}

#' @export
Search_species <- 
function( Class = "predictive", 
          Order = "predictive", 
          Family = "predictive", 
          Genus = "predictive", 
          Species = "predictive",
          add_ancestors = TRUE, 
          Database = FishLife::FishBase_and_RAM, 
          ParentChild_gz = Database$ParentChild_gz ){

  # Match full taxonomy from fishbase
  Match = 1:nrow(rfishbase::fishbase)
  if( Class!="predictive" ) Match = Match[ which(tolower(rfishbase::fishbase$Class[Match])==tolower(Class)) ]
  if( Order!="predictive" ) Match = Match[ which(tolower(rfishbase::fishbase$Order[Match])==tolower(Order)) ]
  if( Family!="predictive" ) Match = Match[ which(tolower(rfishbase::fishbase$Family[Match])==tolower(Family)) ]
  if( Genus!="predictive" ) Match = Match[ which(tolower(rfishbase::fishbase$Genus[Match])==tolower(Genus)) ]
  if( Species!="predictive" ) Match = Match[ which(tolower(rfishbase::fishbase$Species[Match])==tolower(Species)) ]
  if( length(Match)==0 ) stop( paste("Inputs not found in FishBase, please check spelling of",tolower(Class),tolower(Order),tolower(Family),tolower(Genus),tolower(Species)) )

  # add missing taxonomic levels from FishBase if uniquely defined (and throw error if not)
  full_taxonomy = c(Class, Order, Family, Genus, Species)
  if( !all(c(Species)=="predictive") ){
    if( length(unique(rfishbase::fishbase[Match,'Species']))!=1) stop("inputs are not unique")
    if( length(unique(rfishbase::fishbase[Match,'Species']))==1) full_taxonomy[5] = unique(rfishbase::fishbase[Match,'Species'])[1]
  }
  if( !all(c(Species,Genus)=="predictive") ){
    if( length(unique(rfishbase::fishbase[Match,'Genus']))!=1) stop("inputs are not unique")
    if( length(unique(rfishbase::fishbase[Match,'Genus']))==1) full_taxonomy[4] = unique(rfishbase::fishbase[Match,'Genus'])[1]
  }
  if( !all(c(Species,Genus,Family)=="predictive") ){
    if( length(unique(rfishbase::fishbase[Match,'Family']))!=1) stop("inputs are not unique")
    if( length(unique(rfishbase::fishbase[Match,'Family']))==1) full_taxonomy[3] = unique(rfishbase::fishbase[Match,'Family'])[1]
  }
  if( !all(c(Species,Genus,Family,Order)=="predictive") ){
    if( length(unique(rfishbase::fishbase[Match,'Order']))!=1) stop("inputs are not unique")
    if( length(unique(rfishbase::fishbase[Match,'Order']))==1) full_taxonomy[2] = unique(rfishbase::fishbase[Match,'Order'])[1]
  }
  if( !all(c(Species,Genus,Family,Order,Class)=="predictive") ){
    if( length(unique(rfishbase::fishbase[Match,'Class']))!=1) stop("inputs are not unique")
    if( length(unique(rfishbase::fishbase[Match,'Class']))==1) full_taxonomy[1] = unique(rfishbase::fishbase[Match,'Class'])[1]
  }
  match_taxonomy = full_taxonomy

  # Match in database
  Count = 1
  Group = NA
  while( is.na(Group) ){
    Group = match( paste(tolower(match_taxonomy),collapse="_"), tolower(ParentChild_gz[,'ChildName']) )
    if( is.na(Group) ){
      match_taxonomy[length(match_taxonomy)-Count+1] = "predictive"
      Count = Count+1
    }
  }
  message( "Closest match: ", as.character(ParentChild_gz[Group,'ChildName']) )

  # Pick out ancestors
  if( add_ancestors==TRUE ){
    Group = Find_ancestors(child_num=Group, ParentChild_gz=ParentChild_gz)
  }

  # Function to add predictive to taxon name
  Add_predictive = function( char_vec ){
    return_vec = char_vec
    for(i in 1:length(return_vec)){
      vec = strsplit(as.character(return_vec[i]),"_")[[1]]
      return_vec[i] = paste( c(vec,rep("predictive",5-length(vec))), collapse="_")
    }
    return(return_vec)
  }
  match_taxonomy = unique(as.character(Add_predictive(ParentChild_gz[Group,'ChildName'])))
  # Find new matches
  GroupNum = match(match_taxonomy,ParentChild_gz[,'ChildName'])

  # Return match
  Return = list( "GroupNum"=GroupNum, "match_taxonomy"=match_taxonomy )
  return( Return )
}


