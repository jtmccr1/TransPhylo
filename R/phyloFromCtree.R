#' Converts a phylogenetic tree into an ape phylo object with node.comments identifying the host of each node
#' @param ctree combined tree
#' @return phylo object
#' @seealso phyloFromPtree 
#' @export
phyloFromCtree <- function(ctree) {
  nam=ctree$nam
  ntips = length(nam)
  ctree=ctree$ctree
  # each row is a node it can be a tip, an internal bifurcation or an internal node with 1 child representing a host change.
  n<-nrow(ctree)
  if (n==1) return(ape::read.tree(text='(1);'))
  tr<-list()
  tr$Nnode<-n-ntips
  tr$tip.label<-nam
  tr$edge<-matrix(0,n-1,2)
  tr$edge.length<-rep(0,n-1)
  iedge<-1
  
  root<-which(ctree[,1]==min(ctree[,1]))
  tra<-c(1:ntips,root,setdiff((ntips+1):n,root))
  tra0<-tra
  tra2<-1:length(tra)
  tra[tra]<-tra2
  nodeComments = rep(NA,n)
  nodeComments[1:ntips]<-paste0("&host=",nam)
  # tra[i] is the source of each edge. 
  for (i in (ntips+1):n) {
    
    tr$edge[iedge,]<-c(tra[i],tra[ctree[i,2]])
    tr$edge.length[iedge]<-ctree[ctree[i,2],1]-ctree[i,1]
    host<-ctree[i,4]
    
    if(host<=ntips&host>0){
      nodeComments[tra[i]]<-paste0("&host=",nam[host])
    }else{
      nodeComments[tra[i]]<-paste0("&host=","unSampled",host)
    }
    iedge<-iedge+1
    
    if(ctree[i,3]!=0){
      tr$edge[iedge,]<-c(tra[i],tra[ctree[i,3]])
      tr$edge.length[iedge]<-ctree[ctree[i,3],1]-ctree[i,1]
      host<-ctree[i,4]
      
      if(host<=ntips&host>0){
        nodeComments[tra[i]]<-paste0("&host=",nam[host])
      }else{
        nodeComments[tra[i]]<-paste0("&host=","unsampled",host)
      }
      iedge<-iedge+1
    }
  }
  tr$node.comment<-nodeComments
  class(tr)<-'phylo'
  return(tr)
}