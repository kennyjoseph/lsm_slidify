get_block_similarity_matrix <- function(N_ACTORS,N_COVARIATES,N_GROUPS,OUTGROUP_TIE,INGROUP_TIE){
  ##BASELINE_TIE SIMILARITY MATRIX
  similarity_matrix <- matrix(0, nrow=N_ACTORS,ncol=N_ACTORS)
  similarity_matrix[upper.tri(similarity_matrix)] <- runif(N_ACTORS*(N_ACTORS-1)/2,0,OUTGROUP_TIE*2)
  
  groupings <- data.frame(id=1:N_ACTORS)
  for(i in 1:N_COVARIATES){
    groupings[,paste0("Group",i)] <- sample(1:as.numeric(N_GROUPS),
                                            N_ACTORS,replace=T)
  }
  ##Determine matrix of co-memberships in groups, normalized by number of groups
  if(N_COVARIATES > 1){
    percent_shared_memberships <- Reduce("+",
                                         lapply(groupings[,2:ncol(groupings)],function(l){
                                           ifelse(outer(l,l, FUN="-") == 0, 1,0)
                                         }))/(N_COVARIATES)
  } else {
    percent_shared_memberships <- ifelse(outer(groupings[,2],groupings[,2], FUN="-") == 0, 1,0)
  }
  similarity_matrix <- similarity_matrix + INGROUP_TIE*percent_shared_memberships
  diag(similarity_matrix) <- 0; similarity_matrix[lower.tri(similarity_matrix)] <- 0
  similarity_matrix <- ifelse(similarity_matrix > 1, 1, similarity_matrix)
  return(list(similarity_matrix=similarity_matrix, groupings=groupings))
}

get_random_draw_from_matrix <- function(similarity_matrix){
  N_ACTORS <- dim(similarity_matrix)[1]
  random_draw <- matrix(0,nrow=N_ACTORS,ncol=N_ACTORS)
  random_draw[upper.tri(random_draw)] <- rbinom(rep(1,N_ACTORS*(N_ACTORS-1)/2),1,
                                                as.vector(similarity_matrix[upper.tri(similarity_matrix)]))
  lower_indices <- lower.tri(random_draw)
  random_draw[lower_indices] <- t(random_draw)[lower_indices]
  return(random_draw)
}

run_latentnet <- function(net,
                          groupings,
                          n_groups,
                          cov_mult){
  require(latentnet)
  formula <- paste0("net~euclidean(2,G=",n_groups,")")#",var.mul=",cov_mult,")")
  ln_runs <- list(ergmm(as.formula(formula)))
  
  if(N_COVARIATES > 0){
    for(i in 1:N_COVARIATES){
      formula <- paste0(formula,'+nodematch("',colnames(groupings)[i+1],'")')
      z <- ergmm(as.formula(formula))
      ln_runs <- c(ln_runs,list(z))
    }
  }
  return(ln_runs)
}

run_cidnetworks <- function(network,
                            groupings,
                            n_groups,
                            cov_mult){
  if("package:spatstat" %in% search() ){
    detach("package:spatstat", unload=T)
  }
  if("package:CIDnetworks" %in% search()){
    detach("package:CIDnetworks")
  }
  require(CIDnetworks)
  ##Groupings has an ID column
  N_COVARIATES <- ncol(groupings)-1
  
  if(N_COVARIATES > 1){
    edgewise_covariates <- do.call(cbind,lapply(
                                  groupings[,2:ncol(groupings)],function(l){
                                    as.vector(ifelse(outer(l,l, FUN="-") == 0, 1,0))
                                  }))
  } else {
    edgewise_covariates <- as.matrix(as.vector(ifelse(outer(groupings[,2],groupings[,2], FUN="-") == 0, 1,0)),ncol=1)
  }
  ln_runs <- list(CID.Gibbs(sociomatrix=network,components=list(LSM(2))))
  
  ###I can't get this to run with a single covariate-
  ##I get: 
  if(N_COVARIATES > 0){
    for(i in 1:N_COVARIATES){
      if(i == 1){
        e_covariates <- as.matrix(edgewise_covariates,ncol=1)
        ##This is TRUE:
        is.matrix(e_covariates)
        #but from CIDnetworks I get:
        ###Fitting: ordinal outcome with 2 states.
        ###Reinitializing LSM Positions
        ###Error: invalid assignment for reference class field ‘covariates’, should be from class “matrix” or a subclass (was class “numeric”)
        
        ##So lets add a dummy class, just to run something, as its fine w/ 2 columns
        e_covariates <- cbind(e_covariates,1:nrow(e_covariates))
      } else{
        e_covariates <- edgewise_covariates[,1:i]
      }
      z <- CID.Gibbs(sociomatrix=network,components=list(LSM(2),COV(e_covariates)))
      ln_runs <- c(ln_runs,list(z))
    }
  }
  return(ln_runs)
}

get_ripleys_l <- function(ergmm_result){
  require(spatstat)
  d <- ergmm_result$mkl$Z
  p <- as.ppp(d,c(min(d[,1]),max(d[,1]),min(d[,2]),max(d[,2])))
  return(Lest(p, correction="best"))
}


.doFD <- function(data,dimension.type="correlation"){
  dimension.arg = switch(dimension.type,
                         correlation="-q2",
                         hausdorff="-q0")
  FDNQ.tmpfile <- "tmp.dat"
  setwd("/Users/kjoseph/Dropbox/Kenny/papers/current/net_models/fdnq_h/")
  write.table(data,file=FDNQ.tmpfile, row.names=F,col.names=F, sep=" ")
  fdnq.script <- "fdnq.pl"
  cmd = paste(fdnq.script, dimension.arg,   FDNQ.tmpfile)
  str <- system2("perl" ,cmd,stdout=TRUE)
  dat <- read.csv(paste(FDNQ.tmpfile,".points",sep=""), header=FALSE, sep=" ")
  names(dat) = c("mass","bins")
  return(dat)
}

