
#' @export
# Class="predictive"; Order="predictive"; Family="predictive"; Genus="predictive"; Species="predictive"
Search_species = function( Class="predictive", Order="predictive", Family="predictive", Genus="predictive", Species="predictive", add_ancestors=TRUE ){
  # Match full taxonomy from fishbase
  data( fishbase, package="rfishbase")
  Match = 1:nrow(fishbase)
  if( Class!="predictive" ) Match = Match[ which(tolower(fishbase[Match,'Class'])==tolower(Class)) ]
  if( Order!="predictive" ) Match = Match[ which(tolower(fishbase[Match,'Order'])==tolower(Order)) ]
  if( Family!="predictive" ) Match = Match[ which(tolower(fishbase[Match,'Family'])==tolower(Family)) ]
  if( Genus!="predictive" ) Match = Match[ which(tolower(fishbase[Match,'Genus'])==tolower(Genus)) ]
  if( Species!="predictive" ) Match = Match[ which(tolower(fishbase[Match,'Species'])==tolower(Species)) ]
  if( length(Match)==0 ) stop("Inputs not found in FishBase, please check spelling")

  # add missing taxonomic levels from FishBase if uniquely defined (and throw error if not)
  full_taxonomy = c(Class, Order, Family, Genus, Species)
  if( !all(c(Species)=="predictive") ){
    if( length(unique(fishbase[Match,'Species']))!=1) stop("inputs are not unique")
    if( length(unique(fishbase[Match,'Species']))==1) full_taxonomy[5] = unique(fishbase[Match,'Species'])[1]
  }
  if( !all(c(Species,Genus)=="predictive") ){
    if( length(unique(fishbase[Match,'Genus']))!=1) stop("inputs are not unique")
    if( length(unique(fishbase[Match,'Genus']))==1) full_taxonomy[4] = unique(fishbase[Match,'Genus'])[1]
  }
  if( !all(c(Species,Genus,Family)=="predictive") ){
    if( length(unique(fishbase[Match,'Family']))!=1) stop("inputs are not unique")
    if( length(unique(fishbase[Match,'Family']))==1) full_taxonomy[3] = unique(fishbase[Match,'Family'])[1]
  }
  if( !all(c(Species,Genus,Family,Order)=="predictive") ){
    if( length(unique(fishbase[Match,'Order']))!=1) stop("inputs are not unique")
    if( length(unique(fishbase[Match,'Order']))==1) full_taxonomy[2] = unique(fishbase[Match,'Order'])[1]
  }
  if( !all(c(Species,Genus,Family,Order,Class)=="predictive") ){
    if( length(unique(fishbase[Match,'Class']))!=1) stop("inputs are not unique")
    if( length(unique(fishbase[Match,'Class']))==1) full_taxonomy[1] = unique(fishbase[Match,'Class'])[1]
  }
  match_taxonomy = full_taxonomy

  # Match in database
  Count = 1
  Group = NA
  while( is.na(Group) ){
    Group = match( paste(tolower(match_taxonomy),collapse="_"), tolower(ParentChild_gz[,'ChildName']) )
    if( is.na(Group) ) match_taxonomy[length(match_taxonomy)-Count+1] = "predictive"
  }
  message( "Closest match: ", as.character(ParentChild_gz[Group,'ChildName']) )

  # Pick out ancestors
  if( add_ancestors==TRUE ){
    Group = Find_ancestors(Group)
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


