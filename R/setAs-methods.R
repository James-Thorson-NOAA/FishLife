

#' Convert FishLife to other package formats
#'
#' Convert output from package FishLife to package-phylobase or package-sem
#'
#' @name setAs
#' @docType methods
#' @section Usage: \code{as(object, "phylo4d")}
#' @section Usage: \code{as(object, "sem")}
#' @seealso generic \code{\link[methods]{as}}, \code{\link[sem]{sem}} from the
#' \code{sem} package, and \code{\link[phylobase]{phylo4d}} from the
#' \code{phylobase} package.
#' @keywords methods
#' @aliases as as-method as,sem,sem-method
setAs("FishLife", "sem", function(from, to) {

  Sprime = from$Cov_jj
    rownames(Sprime) = colnames(Sprime) = colnames(from$Y_ij)
  out = sem::sem( from$SEM_model,
             S = Sprime,
             N = nrow(from$Y_ij) )   

  # pass out
  return(out)

})

setAs("FishLife", "phylo4d", function(from, to) {

  #
  tip_traits = from$beta_gv[1:ape::Ntip(from$tree),,drop=FALSE]
  rownames(tip_traits) = from$tree$tip.label

  # 
  if(is.null(from$tree$node.label)){
    from$tree$node.label = paste0( "n", 1:ape::Nnode(from$tree) )
  }
  node_traits = from$beta_gv[ape::Ntip(from$tree) + 1:(ape::Nnode(from$tree)-1),,drop=FALSE]
  node_traits = rbind( from$ParHat$alpha_j, node_traits)
  rownames(node_traits) = from$tree$node.label
  out = phylobase::phylo4d( x=from$tree, tip.data=tip_traits, node.data=node_traits )

  # pass out
  return(out)

})

