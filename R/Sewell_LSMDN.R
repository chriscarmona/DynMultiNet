#' @title
#'    Bayesian Learning of Dynamic Networks
#'
#' @description
#'    \code{LSMDN} Implements model from Sewell and Chen, 2015
#'
#' @param Y Array. Network information.
#' @param p Integer.  Dimension of the Euclidean latent space
#' @param MissData Boolean. Are there missing edges?
#' @param N Integer. Number of total MCMC iterations.
#' @param burnin Integer. Number of MCMC iterations in the warm-up period.
#' @param n0 Integer. Size of subsample when llApprox==TRUE.
#' @param llApprox Boolean. Use log likelihood approximation?
#' @param tuneX real. MCMC tuning parameters
#' @param tuneBIO real. MCMC tuning parameters
#' @param Kappa real. MCMC tuning parameters
#' @details
#'    The model assumes a latent variable approach
#'
#' @return
#'    A list with the following components:
#' \describe{
#'     \item{\code{theta_mcmc}}{Matrix with the chain of the parameters in the model.}
#' }
#'
#' @useDynLib DynMultiNet, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' 
#' @import MCMCpack
#' @importFrom vegan procrustes
#' 
#' @export
#' 

LSMDN <- function( Y, p=2,
                   MissData=FALSE,
                   N=20000, burnin = round(N/10), n0=100,
                   llApprox=FALSE,
                   tuneX=0.0075, tuneBIO=0.1, Kappa=175000,
                   out_rds=NULL
) {
  
  if(!llApprox) n0<-NULL
  
  # This software only deal with binary edges
  Y_orig <- Y
  Y[Y>0] <- 1
  
  
  
  #If MissData==TRUE, construct Missing: Missing[[tt]]= c(1,2,3) => at time tt we have no row data on actors 1,2&3 
  if(MissData){
    Missing <- list() #Enter as list manually, or if NAs are in data run the following:
    temp = which(is.na(Y),arr.ind=TRUE)
    for(tt in 1:dim(Y)[3]){
      Missing[[tt]] = unique(temp[which(temp[,3]==tt),1]) # which rows have missing data for each time
    }
    Y[temp] = 0;rm(temp)
  }
  
  pb <- txtProgressBar(min=2,max=N,style=3)
  on.exit(close(pb))
  
  #######################
  ###     Startup     ###
  #######################
  
  n <- dim(Y)[1]
  TT <- dim(Y)[3]
  
  ###
  ###Missing Data
  ###
  if(MissData){
    Yaggreg <- array(0,dim(Y[,,1]))
    for(tt in 1:TT) Yaggreg <- Yaggreg + Y[,,tt]
    Yaggreg[which(Yaggreg>1,arr.ind=TRUE)] <- 1
    SPaggreg <- shortest.paths(graph=graph.adjacency(Yaggreg),mode="all")
    SPaggreg[which(SPaggreg==Inf,arr.ind=TRUE)] <- 5
    
    outDeg <- inDeg <- denom <- numeric(n)
    for(tt in 1:TT){
      denom[c(1:n)[-Missing[[tt]]]] <- denom[c(1:n)[-Missing[[tt]]]] + 1
      if(length(Missing[[tt]])==0) denom <- denom + 1
      inDeg <- inDeg + colSums(Y[,,tt])/(n-1-length(Missing[[tt]])+as.numeric(1:n%in%Missing[[tt]]))*(n-1)
      outDeg <- outDeg + rowSums(Y[,,tt])
    }
    inDeg= inDeg/TT
    outDeg <- round(outDeg/denom)
    
    for(tt in 1:TT){
      for(i in Missing[[tt]]){
        Probs <- inDeg[-i]/SPaggreg[i,-i]
        # ind <- sample(size=outDeg[i],x=c(1:n)[-i],prob=Probs,replace=FALSE) # this doesn't work
        ind <- sample(size=rowSums(Y[,,tt]),x=c(1:n)[-i],prob=Probs,replace=FALSE)
        Y[ind,i,tt] <- Y[i,ind,tt]  <- 1
      }
    }
  }
  
  Bout <- Bin <- numeric(N)
  w <- matrix(0,n,N)
  s2 <- t2 <- numeric(N)
  # AccRate <- matrix(0,N,3+n*TT); AccRate[1,] <- rep(NA,3+n*TT)
  # colnames(AccRate) <- c("Bin","Bout","weights",
  #                        paste("X",rep(1:n,TT),rep(1:TT,each=n),sep=","))
  AccRate <- numeric(3+n*TT)
  names(AccRate) <- c("Bin","Bout","weights",
                      paste("X",rep(1:n,TT),rep(1:TT,each=n),sep=","))
  
  
  ###
  ###Weights
  ###
  for(tt in 1:TT){
    w[,1] <- w[,1] + apply(Y[,,tt],1,sum) + 
      apply(Y[,,tt],2,sum)
  }
  w[,1] <- w[,1]/sum(Y)/2
  if(sum(w==0)>0){
    w[,1] <- w[,1]+1e-5
    w[,1] <- w[,1]/sum(w[,1])
  }
  w[,1] <- w[,1]/sum(w[,1])
  
  ###
  ###Initial Latent Positions (GMDS, Sarkar and Moore, 2005)
  ###
  dissim <- array(0,dim=dim(Y))
  for(tt in 1:TT) dissim[,,tt] <- shortest.paths(graph=graph.adjacency(Y[,,tt]),mode="all")
  dissim[which(dissim==Inf,arr.ind=TRUE)] <- max(c(dissim[which(dissim!=Inf,arr.ind=TRUE)]))
  X <- list()
  X[[1]] <- array(0,c(p,TT,n))
  X[[1]][,1,] <- t(cmdscale(d=dissim[,,1],k=p))
  temp.lambda <- 10
  H <- matrix(-1/n,n,n)+diag(n)
  
  for(tt in 2:TT){
    temp <- 1/(1+temp.lambda)*H%*%(-1/2*dissim[,,tt]^2)%*%H +
      temp.lambda/(1+temp.lambda)*t(X[[1]][,tt-1,])%*%X[[1]][,tt-1,]
    temp <- eigen(temp)
    X[[1]][,tt,] <- t(temp$vectors[,1:p]%*%diag(temp$values[1:p])^(1/2))
    X[[1]][,tt,] <- t(vegan::procrustes(X=t(X[[1]][,tt-1,]),Y=t(X[[1]][,tt,]),scale=FALSE)$Yrot)
  }
  
  ###
  ###Initial \beta_{IN} and \beta_{OUT}; scale latent positions
  ###
  initialize.wrap <- function(x){
    -c_initialize1(X[[1]],c(n,p,TT),Y,XSCALE=1/n,
                   BETAIN=x[1],BETAOUT=x[2],w[,1])
  }
  initialize.grad.wrap <- function(x){
    -c_initialize1_grad(X[[1]],c(n,p,TT),Y,XSCALE=1/n,
                        BETAIN=x[1],BETAOUT=x[2],w[,1])[2:3]
  }
  Optim <- optim(par=c(1,1),fn=initialize.wrap,
                 gr=initialize.grad.wrap,method="BFGS")
  X[[1]] <- X[[1]]/n
  Bin[1] <- max(Optim$par[1],1e-4)
  Bout[1] <- max(Optim$par[2],1e-4)
  
  Xit0 <- t(X[[1]][,1,])
  for(tt in 2:TT)Xit0 <- rbind(Xit0,t(X[[1]][,tt,]))
  Xit0 <- Xit0 -
    kronecker(rep(1,n*TT),matrix(apply(Xit0,2,mean),1,p))
  
  ###
  ###Priors and further initializations:
  ###
  nuIn <- Bin[1]
  nuOut <- Bout[1]
  xiIn <- xiOut <- 100
  
  #Tau^2
  t2[1] <- sum(X[[1]][,1,]*X[[1]][,1,])/(n*p)
  shapeT2 <- 2.05
  scaleT2 <- (shapeT2-1)*t2[1]
  
  #Sigma^2
  s2[1] <- 0.001
  shapeS2 <- 9
  scaleS2 <- 1.5
  
  
  ###
  ###log-likelihood approximation subsampling
  ###
  if(llApprox){
    DEGREE <- array(0,c(n,2,TT))
    for(tt in 1:TT){
      DEGREE[,1,tt] <- colSums(Y[,,tt])
      DEGREE[,2,tt] <- rowSums(Y[,,tt])
    }
    dinmax <- max(DEGREE[,1,])
    doutmax <- max(DEGREE[,2,])
    ELIN <- array(0,c(n,dinmax,TT))
    ELOUT <- array(0,c(n,doutmax,TT))
    
    for(tt in 1:TT){
      for(i in 1:n){
        ind <- which(Y[,i,tt]==1)
        if(length(ind)>0){
          ELIN[i,1:length(ind),tt] <- ind
        }
        ind <- which(Y[i,,tt]==1)
        if(length(ind)>0){
          ELOUT[i,1:length(ind),tt] <- ind
        }
      }
    }
    
    n0<-max(n0,dinmax+10,doutmax+10)
    
    edgeList <- list()
    for(i in 1:n){
      edgeList[[i]] <- which(Y[i,,]==1|Y[,i,]==1,arr.ind=TRUE)[,1]
      edgeList[[i]] <- unique(edgeList[[i]])
    }
    
    SUBSEQ = matrix(0,n,n0)
    for(i in 1:n){
      #   SUBSEQ[i,] <- sample(c(1:n)[-i],size=n0,replace=FALSE) #Simple random sample
      nOnes <- round(length(edgeList[[i]])/n*n0) #stratified sampling
      if(length(edgeList[[i]])>0){ nOnes <- max(nOnes,1) }
      SUBSEQ[i,1:nOnes] <- sample(edgeList[[i]],size=nOnes,replace=FALSE)
      SUBSEQ[i,(nOnes+1):n0] <- sample(c(1:n)[-c(i,edgeList[[i]])],size=n0-nOnes,replace=FALSE)
    }
  }
  
  
  ##############################
  ###     LSMDN SAMPLING     ###
  ##############################
  
  for(it in 2:N){ # it <- 2
    #   for(it in 2:500){
    
    RN <- rnorm(n*TT*p)
    RNBIO <- rnorm(2)
    if(llApprox){
      if(it%%100==0){
        SUBSEQ = matrix(0,n,n0)
        for(i in 1:n){
          nOnes <- round(length(edgeList[[i]])/n*n0) #stratified sampling
          if(length(edgeList[[i]])>0){ nOnes <- max(nOnes,1) }
          SUBSEQ[i,1:nOnes] <- sample(edgeList[[i]],size=nOnes,replace=FALSE)
          SUBSEQ[i,(nOnes+1):n0] <- sample(c(1:n)[-c(i,edgeList[[i]])],size=n0-nOnes,replace=FALSE)
        }
      }
    }
    if(llApprox){
      Draws <- c_update2(X[[it-1]],c(n,p,TT,dinmax,doutmax),tuneX,Y,
                         Bin[it-1],Bout[it-1],tuneBIO,w[,it-1],
                         t2[it-1],s2[it-1],xiIn,xiOut,nuIn,
                         nuOut,CAUCHY=0,RN,RNBIO,ELOUT,ELIN,SUBSEQ,DEGREE)
    }else{
      Draws <- c_update1(X[[it-1]],c(n,p,TT,1),tuneX,Y,
                         Bin[it-1],Bout[it-1],tuneBIO,w[,it-1],
                         t2[it-1],s2[it-1],xiIn,xiOut,nuIn,
                         nuOut,CAUCHY=0,RN,RNBIO)
    }
    X[[it]] <- Draws[[1]]
    Bin[it] <- Draws[[2]]
    Bout[it] <- Draws[[3]]
    AccRate <- AccRate + Draws[[4]]
    
    if(it==burnin){
      Xit0 <- t(X[[it]][,1,]) # X[[it]] is an array with dim (dim latent,time,num actors)
      for(tt in 2:TT) Xit0 <- rbind(Xit0,t(X[[it]][,tt,]))
    }
    if(it>burnin){
      XitCentered <- t(X[[it]][,1,])
      for(tt in 2:TT) XitCentered <- rbind(XitCentered,t(X[[it]][,tt,]))
      procr <- vegan::procrustes(X=Xit0,Y=XitCentered,scale=FALSE)$Yrot
      for(tt in 1:TT){
        X[[it]][,tt,] <- t(procr[((tt-1)*n+1):(tt*n),])
      }
    }
    if(it < N) X[[it+1]] <- X[[it]]
    
    #------------------Step 2--------------------------------
    Draws1 <- c_t2s2Parms(X[[it]],c(n,p,TT,1),shapeT2,
                          shapeS2,scaleT2,scaleS2)
    t2[it] <- rinvgamma(1,shape=Draws1[[1]],scale=Draws1[[2]])
    s2[it] <- rinvgamma(1,shape=Draws1[[3]],scale=Draws1[[4]])
    
    #------------------Step 3--------------------------------
    
    w[,it] <- rdirichlet(1,alpha=Kappa*w[,it-1])
    if(llApprox){
      Draws2 <- c_WAccProb2(X[[it]],c(n,p,TT,dinmax,doutmax),Y,
                            Bin[it],Bout[it],Kappa,w[,it-1],w[,it],
                            ELOUT,ELIN,SUBSEQ,DEGREE)
    }else{
      Draws2 <- c_WAccProb1(X[[it]],c(n,p,TT,1),Y,
                            Bin[it],Bout[it],Kappa,w[,it-1],w[,it])
    }
    w[,it] <- Draws2[[1]]
    AccRate[3] <- AccRate[3] + Draws2[[2]]
    
    #------------------Step 4--------------------------------
    
    if(MissData){
      for(tt in 1:TT){
        Y <- c_missing(X[[it]],c(n,p,TT),MMM=Missing[[tt]]-1,Y,Ttt=tt,
                       BETAIN=Bin[it],BETAOUT=Bout[it],WW=w[,it])
      }
    }
    
    
    if( !is.null(out_rds) & is.element(it, floor(N*seq(0,1,0.25)[-1]) ) ) {
      out <- list( Y=Y_orig,
                   p=p, TT=TT, n=n,
                   MissData=MissData,
                   X=X, # Latent coordinates
                   N=N, burnin=burnin, n0=n0,
                   Bin=Bin, Bout=Bout, # beta_in and beta_out
                   s2=s2, # sigma
                   t2=t2, # tau
                   w=w, # weights
                   AccRate=AccRate )
      saveRDS(out,file=out_rds)
    }
    
    setTxtProgressBar(pb,it)
  }
  
  out <- list( Y=Y_orig,
               p=p, TT=TT, n=n,
               MissData=MissData,
               X=X, # Latent coordinates
               N=N, burnin=burnin, n0=n0,
               Bin=Bin, Bout=Bout, # beta_in and beta_out
               s2=s2, # sigma
               t2=t2, # tau
               w=w, # weights
               AccRate=AccRate )
  
  return(out)
  
}

#' @export
c_postzeroprob <- function(Xi1, Xi2, Xj1, Xj2, SS2, LAM, PP0) {
  c_postzeroprob(Xi1, Xi2, Xj1, Xj2, SS2, LAM, PP0)
}

#' @export
c_prediction <- function(EX, SIG2, X1T, X2T, BIN, BOUT, WW) {
  c_prediction(EX, SIG2, X1T, X2T, BIN, BOUT, WW)
}