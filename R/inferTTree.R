#' A helper function to start the tree log of a run
#' It writes the taxa labels from a provided phylo object, but does not add any trees.
#' 
#' @param file the log file
#' @param tree phylo object which provides taxa labels and encodings
startNexus<-function(file,tree){
  string = "#NEXUS\n\nBegin taxa;\n\tDimensions ntax="
  string = paste0(string,length(tree$tip.label),";\n\tTaxlabels\n")
  for(i in 1:length(tree$tip.label)){
    string = paste0(string,"\t\t",tree$tip.label[i],"\n")
  }
  string=paste0(string,"\t\t;\nEND;\n\nBegin trees;\n\n\tTranslate\n")
  for(i in 1:(length(tree$tip.label)-1)){
    string = paste0(string,"\t\t",i," ",tree$tip.label[i],",\n")
  }
  lastTaxa<-length(tree$tip.label)
  string = paste0(string,"\t\t",lastTaxa," ",tree$tip.label[lastTaxa],"\n\t\t;")
  string=paste0(string,"\nEND;\n")
  cat(string,file=file)
}


#' This appends a tree to a nexus file
#' This assumes the taxa labels are already in the file. 
#' The last line of the file is assumed to be "END;". It is removed 
#' and "END;" is added after the tree.
#' 
#' @param file The tree file we are appending to
#' @param tree The phylo tree object
#' @param name The name of the tree file.
appendTree<-function(file,tree,name="unnamed"){
  fileString = readLines(file)
  fileString=fileString[-length(fileString)] # get rid of the last line "END;"
  if(is.rooted(tree)){
    rootedLab = "[&R] "
  }else{
    rootedLab = "[&U] "
  }
  treeString<-paste0("tree ",name, " = ",rootedLab, print_annotated(tree,"newick"),";\nEND;")
  fileString[length(fileString)+1]=treeString
  writeLines(fileString,file)
}
#' A helper function to start the log of a run
#' It just writes the header
#' @param file the log file
startLog<-function(file){
  logString = paste(c("state","pTTree","pPTree","neg","off.r","off.p","pi","w.shape","w.scale","ws.shape","ws.scale\n"),collapse="\t")
  cat(logString,file=file)
}
#'Append States to a log file
#' @param file the file to append to
#' @param record named list that holds the parameters we are logging
#' @param state the integer number of the state
appendLog<-function(file,record,state){
  logItems = c(state,unlist(record[!(names(record) %in% c("ctree","source"))]))
  logString = paste(logItems,collapse = "\t")
  logString=paste0(logString,"\n")
  cat(logString,file=file,append=T)
}




#' Infer transmission tree given a phylogenetic tree
#' @param ptree Phylogenetic tree
#' @param fileRoot The file root to log states.Defualt is NULL - no logging
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation time
#' @param w.scale Scale parameter of the Gamma probability density function representing the generation time 
#' @param ws.shape Shape parameter of the Gamma probability density function representing the sampling time
#' @param ws.scale Scale parameter of the Gamma probability density function representing the sampling time 
#' @param mcmcIterations Number of MCMC iterations to run the algorithm for
#' @param thinning MCMC thinning interval between two sampled iterations
#' @param startNeg Starting value of within-host coalescent parameter Ne*g
#' @param startOff.r Starting value of parameter off.r
#' @param startOff.p Starting value of parameter off.p
#' @param startPi Starting value of sampling proportion pi
#' @param updateNeg Whether of not to update the parameter Ne*g
#' @param updateOff.r Whether or not to update the parameter off.r
#' @param updateOff.p Whether or not to update the parameter off.p
#' @param updatePi Whether or not to update the parameter pi
#' @param startCTree Optional combined tree to start from
#' @param updateTTree Whether or not to update the transmission tree
#' @param optiStart Whether or not to optimise the MCMC start point
#' @param dateT Date when process stops (this can be Inf for fully simulated outbreaks)
#' @return posterior sample set of transmission trees
#' @import phylotate
#' @import ape
#' @export
inferTTree = function(ptree, fileRoot=NULL, w.shape=2, w.scale=1, ws.shape=w.shape, ws.scale=w.scale, mcmcIterations=1000,
                      thinning=1, startNeg=100/365, startOff.r=1, startOff.p=0.5, startPi=0.5, updateNeg=TRUE,
                      updateOff.r=TRUE, updateOff.p=FALSE, updatePi=TRUE, startCTree=NA, updateTTree=TRUE,
                      optiStart=TRUE, dateT=Inf) {
#  memoise::forget(getOmegabar)
#  memoise::forget(probSubtree)
  if(!(is.null(fileRoot))){
    ## Set up loggers
    logFile = paste0(fileRoot,".log")
    # treesfile = paste0(fileRoot,".trees")
    #Start logs
    startLog(logFile)
    # phyloTree<-phyloFromPTree(ptree) # just to get the tip names ect.
    # startNexus(treesfile,phyloTree)
  }
  ptree$ptree[,1]=ptree$ptree[,1]+runif(nrow(ptree$ptree))*1e-10#Ensure that all leaves have unique times
  for (i in (ceiling(nrow(ptree$ptree)/2)+1):nrow(ptree$ptree)) for (j in 2:3) 
    if (ptree$ptree[ptree$ptree[i,j],1]-ptree$ptree[i,1]<0) 
      stop("The phylogenetic tree contains negative branch lengths!")
  
  #MCMC algorithm
  neg <- startNeg
  off.r <- startOff.r
  off.p <- startOff.p
  pi <- startPi
  if (is.na(sum(startCTree))) ctree <- makeCtreeFromPTree(ptree,ifelse(optiStart,off.r,NA),off.p,neg,pi,w.shape,w.scale,ws.shape,ws.scale,dateT)#Starting point 
  else ctree<-startCTree
  ttree <- extractTTree(ctree)
  record <- vector('list',mcmcIterations/thinning)
  pTTree <- probTTree(ttree$ttree,off.r,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,dateT) 
  pPTree <- probPTreeGivenTTree(ctree,neg) 
  pb <- txtProgressBar(min=0,max=mcmcIterations,style = 3)
  for (i in 1:mcmcIterations) {#Main MCMC loop
    if (i%%thinning == 0) {
      #Record things 
      setTxtProgressBar(pb, i)
      #message(sprintf('it=%d,neg=%f,off.r=%f,off.p=%f,pi=%f,Prior=%e,Likelihood=%f,n=%d',i,neg,off.r,off.p,pi,pTTree,pPTree,nrow(extractTTree(ctree))))
      record[[i/thinning]]$ctree <- ctree
      record[[i/thinning]]$pTTree <- pTTree 
      record[[i/thinning]]$pPTree <- pPTree 
      record[[i/thinning]]$neg <- neg 
      record[[i/thinning]]$off.r <- off.r
      record[[i/thinning]]$off.p <- off.p
      record[[i/thinning]]$pi <- pi
      record[[i/thinning]]$w.shape <- w.shape
      record[[i/thinning]]$w.scale <- w.scale
      record[[i/thinning]]$ws.shape <- ws.shape
      record[[i/thinning]]$ws.scale <- ws.scale
      record[[i/thinning]]$source <- ctree$ctree[ctree$ctree[which(ctree$ctree[,4]==0),2],4]
      if (record[[i/thinning]]$source<=length(ctree$nam)) record[[i/thinning]]$source=ctree$nam[record[[i/thinning]]$source] else record[[i/thinning]]$source='Unsampled'
    
      if(!(is.null(fileRoot))){
        #Write to logs
        # PhyloCtree<-phyloFromCtree(ctree)
        # appendTree(treesfile,PhyloCtree,paste0("STATE_",i))
        appendLog(logFile,record[[i/thinning]],i)
      }
      
      }
    
    if (updateTTree) {
    #Metropolis update for transmission tree 
    prop <- proposal(ctree$ctree) 
    ctree2 <- list(ctree=prop$tree,nam=ctree$nam)
    ttree2 <- extractTTree(ctree2)
    pTTree2 <- probTTree(ttree2$ttree,off.r,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,dateT) 
    pPTree2 <- probPTreeGivenTTree(ctree2,neg) 
    if (log(runif(1)) < log(prop$qr)+pTTree2 + pPTree2-pTTree-pPTree)  { 
      ctree <- ctree2 
      ttree <- ttree2
      pTTree <- pTTree2 
      pPTree <- pPTree2 
    } 
    }
    
    if (updateNeg) {
      #Metropolis update for Ne*g, assuming Exp(1) prior 
      neg2 <- abs(neg + (runif(1)-0.5)*0.5)
      pPTree2 <- probPTreeGivenTTree(ctree,neg2) 
      if (log(runif(1)) < pPTree2-pPTree-neg2+neg)  {neg <- neg2;pPTree <- pPTree2} 
    }
    
    if (updateOff.r) {
      #Metropolis update for off.r, assuming Exp(1) prior 
      off.r2 <- abs(off.r + (runif(1)-0.5)*0.5)
      pTTree2 <- probTTree(ttree$ttree,off.r2,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,dateT) 
      if (log(runif(1)) < pTTree2-pTTree-off.r2+off.r)  {off.r <- off.r2;pTTree <- pTTree2}
    }
    
    if (updateOff.p) {
      #Metropolis update for off.p, assuming Unif(0,1) prior 
      off.p2 <- abs(off.p + (runif(1)-0.5)*0.1)
      if (off.p2>1) off.p2=2-off.p2
      pTTree2 <- probTTree(ttree$ttree,off.r,off.p2,pi,w.shape,w.scale,ws.shape,ws.scale,dateT) 
      if (log(runif(1)) < pTTree2-pTTree)  {off.p <- off.p2;pTTree <- pTTree2}
    }

        if (updatePi) {
      #Metropolis update for pi, assuming Unif(0.01,1) prior 
      pi2 <- pi + (runif(1)-0.5)*0.1
      if (pi2<0.01) pi2=0.02-pi2
      if (pi2>1) pi2=2-pi2
      pTTree2 <- probTTree(ttree$ttree,off.r,off.p,pi2,w.shape,w.scale,ws.shape,ws.scale,dateT) 
      if (log(runif(1)) < pTTree2-pTTree)  {pi <- pi2;pTTree <- pTTree2}       
    }
    
  }#End of main MCMC loop
  
  #close(pb)
  return(record)
}
