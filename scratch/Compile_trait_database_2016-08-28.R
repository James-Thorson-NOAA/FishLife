
##################
#
# Instructions:
#  1.  Use R v3.3.0 "Supposedly educational" (not revolutionR version)
#
##################


if( !"rfishbase" %in% installed.packages()) install.packages("rfishbase")
library(rfishbase)

##################
# Example: Scrape rockfishes
##################

# Identify rockfishes
species_names = fishbase[ which(fishbase$Family=="Sebastidae"), c("Family","Genus","Species") ]    # apply( species_names, MARGIN=1, FUN=function(vec){paste(rev(rev(vec)[1:2]),collapse=" ")})
species_names = data.frame(species_names, "sciname"=apply( species_names, MARGIN=1, FUN=function(vec){paste(rev(rev(vec)[1:2]),collapse=" ")}) )

# growth curve / mortality
Growth = popgrowth( as.character(species_names[,'sciname']), fields=c("sciname","Sex","Loo","TLinfinity","K","to","Winfinity","tmax","tm","M","Lm") )

##################
# Scrape all fishes
##################

RootFile = "C:/Users/James.Thorson/Desktop/UW Hideaway/Collaborations/2016 -- fish traits on taxonomy/"
DateFile = paste0( RootFile, "Database_", "2018-06-21", "/" )
  dir.create(DateFile)

Family_set = unique(fishbase$Family)

# Loop through database
for( sI in 1:length(Family_set)){
  if( !file.exists(paste0(DateFile,Family_set[sI],".RData")) ){
    # Identify species
    species_names = fishbase[ which(fishbase$Family==Family_set[sI]), c("Class","Order","Family","Genus","Species") ]    # apply( species_names, MARGIN=1, FUN=function(vec){paste(rev(rev(vec)[1:2]),collapse=" ")})
    species_names = data.frame(species_names, "sciname"=apply( species_names, MARGIN=1, FUN=function(vec){paste(rev(rev(vec)[1:2]),collapse=" ")}) )

    # growth curve / mortality
    Fields2Get = c("sciname","Sex","Loo","TLinfinity","K","to","Winfinity","tmax","tm","M","Lm","Type","Temperature")
    Growth = popgrowth( as.character(species_names[,'sciname']), fields=Fields2Get )
    if( nrow(Growth)==0) Growth = matrix(NA, nrow=0, ncol=length(Fields2Get), dimnames=list(NULL,Fields2Get))

    # Match to taxonomy
    DF = cbind( species_names[match(Growth[,'sciname'],species_names[,'sciname']),c("Class","Order","Family","Genus","Species")], Growth )
    save( DF, file=paste0(DateFile,Family_set[sI],".RData"))
  }
}

##################
# Combine data
##################

Family_set = list.files( DateFile )
  #Family_set = sapply(Family_set, FUN=function(char){gsub(x=char, pattern=".RData", replacement="")})

# Loop through database
CombinedDF = NULL
for( sI in 1:length(Family_set)){
  load( file=paste0(DateFile,Family_set[sI]) )
  CombinedDF = rbind( CombinedDF, DF )
}

# Removing non-TL length measurements
CombinedDF[ which(CombinedDF[,'Type']!="TL"), c("Loo","Lm") ] = NA

# Trimming columns
CombinedDF = CombinedDF[,c("Class","Order","Family","Genus","Species","Loo","K","Winfinity","tmax","tm","M","Lm","Temperature")]

# Adding predictive levels for phylogeny
#CombinedDF = rbind( matrix(NA,nrow=1,ncol=ncol(CombinedDF),dimnames=list(NULL,colnames(CombinedDF))), CombinedDF )
#CombinedDF[1,1:5] = "Predictive"

# Save
save( CombinedDF, file=paste0(RootFile,"CombinedDF.RData"))

# Explore
apply( CombinedDF, MARGIN=2, FUN=function(vec){length(unique(vec))})

