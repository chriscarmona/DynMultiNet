#' @title
#'    Synthetic network data
#'
#' @description
#'    \code{gen_synth_net_icml} Generates a multilayer network
#'    
#' @param directed Boolean. Indicates if the provided network is directed, i.e. the adjacency matrix is assymetrical.
#' @param weighted Boolean. Indicates if the provided network is weighted, i.e. edges with values other that 0 and 1.
#' @param rds_file String. Indicates a file (.rds) where the output will be saved.
#' 
#' @export
gen_synth_net_icml <- function( directed=TRUE,
                                weighted=TRUE,
                                rds_file=NULL ) {
  
  # Network topology
  V_net<-15 # number of agents
  T_net<-30 # number of timesteps
  K_net<-2 # number of layers
  P_pred <- 2 # number of predictors
  times_net <- sort(runif(n=T_net,0,T_net)) # times of observation
  
  nodes_net=1:V_net
  layers_net=1:K_net
  
  sigma_k <- c(0.5,0.3)
  
  # Trajectories for the systemic changes in connectivity
  baseline_1<-function(t,t0,t_smooth,a,b1,b2) {a+b1*exp(-((t-t0)/t_smooth)^2)+b2*plogis((t-t0)/t_smooth)}
  baseline_2<-function(t,t0,t_smooth,a,b) {a+b*plogis((t-t0)/t_smooth)}
  
  # Probability of connection
  lambda_ijtk <- array(NA,dim=c(V_net,V_net,T_net,K_net))
  
  # layer 1 #
  # three groups, high connectivity between central agents, low in the perophery #
  lambda_ijtk[,,,1] <- 0.05
  lambda_ijtk[5:8,,,1] <- lambda_ijtk[5:8,,,1] + 0.15 # agents 5:8 are sociable
  lambda_ijtk[,5:8,,1] <- lambda_ijtk[,5:8,,1] + 0.20 # agents 5:8 are popular
  lambda_ijtk[5:8,5:8,,1] <- lambda_ijtk[5:8,5:8,,1] + 0.05 # agents 5:8 connect between them
  lambda_ijtk[1:4,,,1] <- lambda_ijtk[1:4,,,1] + 0.3 # agents 1:4 are more sociable
  lambda_ijtk[,1:4,,1] <- lambda_ijtk[,1:4,,1] + 0.4 # agents 1:4 are more popular
  lambda_ijtk[1:4,1:4,,1] <- lambda_ijtk[1:4,1:4,,1] + 0.15 # agents 1:4 connect more between them
  # lambda_ijtk[,,1,1]
  
  # layer 2 #
  # two groups, more balanced #
  lambda_ijtk[,,,2] <- 0.2
  lambda_ijtk[1:5,,,2] <- lambda_ijtk[1:5,,,2] + 0.2 # agents 1:5 are sociable
  lambda_ijtk[,1:5,,2] <- lambda_ijtk[,1:5,,2] + 0.3 # agents 1:5 are popular
  lambda_ijtk[1:5,1:5,,2] <- lambda_ijtk[1:5,1:5,,2] + 0.2 # agents 1:5 connect more between them
  # lambda_ijtk[,,1,2]
  lambda_ijtk[] <- c(lambda_ijtk)+runif(length(c(lambda_ijtk)),-0.1,0.1)
  lambda_ijtk[] <- pmax(0.01,c(lambda_ijtk)); lambda_ijtk[] <- pmin(0.99,c(lambda_ijtk))
  
  # Linear predictors #
  gamma_ijtk <- lambda_ijtk
  gamma_ijtk[] <- qlogis(c(lambda_ijtk))
  # lambda_ijtk[] <- plogis(gamma_ijtk)
  
  mu_ijtk <- gamma_ijtk
  for(k in 1:K_net) {
    mu_ijtk[,,,k] <- scale(c(mu_ijtk[,,,k])) + c(8,6)[k]
  }
  # hist(mu_ijtk[,,,2],20)
  
  # Baseline theta_t for each layer #
  # curve( baseline_1(x,0.3*T_net,0.2*T_net,0.5,1.5,-1.5),from = 0,to = T_net )
  # curve( baseline_2(x,0.5*T_net,0.2*T_net,0.5,-1.5),from = 0,to = T_net,add=T)
  theta_t <- cbind( baseline_1(times_net,0.3*T_net,0.2*T_net,0.5,1.5,-1.5),
                    baseline_2(times_net,0.5*T_net,0.2*T_net,0.5,-1.5) )
  
  # Adding random dynamics for each node
  # Prior Covariance matrix for GPs
  GP_cov_prior <- outer( times_net, times_net, FUN=function(x,y,k=0.1,delta=0.3*T_net){ k*exp(-((x-y)/delta)^2) } )
  diag(GP_cov_prior) <- diag(GP_cov_prior)+0.001
  # cov2cor(GP_cov_prior) # underlying correlations
  gp_ti_sp_link <- t(mvtnorm::rmvnorm( n = V_net*K_net,
                                       mean = rep(0,T_net),
                                       sigma = GP_cov_prior ))
  # range(gp_ti_sp_link)
  for( k in 1:K_net ) { #k<-1
    for( i in 1:V_net ) { #i<-1
      gamma_ijtk[i,,,k] <- t(t(gamma_ijtk[i,,,k])+gp_ti_sp_link[,(k-1)*V_net+i]+theta_t[,k])
      gamma_ijtk[,i,,k] <- t(t(gamma_ijtk[,i,,k])+gp_ti_sp_link[,(k-1)*V_net+i]+theta_t[,k])
    }
  }
  # new probabilities
  lambda_ijtk[] <- plogis(gamma_ijtk)
  
  # weights #
  # Adding random dynamics for each node #
  GP_cov_prior <- outer( times_net, times_net, FUN=function(x,y,k=0.3,delta=0.3*T_net){ k*exp(-((x-y)/delta)^2) } )
  diag(GP_cov_prior) <- diag(GP_cov_prior)+0.001
  # cov2cor(GP_cov_prior) # underlying correlations
  gp_ti_sp_weight <- t(mvtnorm::rmvnorm( n = V_net*K_net,
                                         mean = rep(0,T_net),
                                         sigma = GP_cov_prior ) )
  # range(gp_ti_sp_weight)
  eta_i <- cbind( baseline_2(times_net,0.3*T_net,0.2*T_net,-1,3),
                  baseline_2(times_net,0.5*T_net,0.2*T_net,1,-3) )
  for( k in 1:K_net ) { #k<-1
    for( i in 1:V_net ) { #i<-1
      # Baseline eta_t for node kind #
      mu_ijtk[i,,,k] <- t(t(mu_ijtk[i,,,k])+gp_ti_sp_weight[,(k-1)*V_net+i]+eta_i[,ifelse( is.element(i,1:5),1,2)])
      mu_ijtk[,i,,k] <- t(t(mu_ijtk[,i,,k])+gp_ti_sp_weight[,(k-1)*V_net+i]+eta_i[,ifelse( is.element(i,1:5),1,2)])
    }
  }
  # mu_ijtk[,,1,1]
  
  # i<-1;j<-11;k<-1
  # plot(y=mu_ijtk[i,j,,k],x=times_net,main=paste("mu i=",i,", j=",j,", k=",k,sep=""))
  # plot(y=gamma_ijtk[i,j,,k],x=times_net,main=paste("gamma i=",i,", j=",j,", k=",k,sep=""))
  
  for( k in 1:K_net ) { #k<-1
    for( t in 1:T_net ) { #i<-1
      diag(lambda_ijtk[,,t,k]) <- diag(mu_ijtk[,,t,k]) <- 0
    }
  }
  
  ### NETWORK ###
  y_ijtk <- array( data=NA, dim=c(V_net,V_net,T_net,K_net) )
  dimnames(y_ijtk) <- list(nodes_net,nodes_net,times_net,layers_net)
  # indicator y_ijtk>0
  z_ijtk <- y_ijtk
  
  ### link incidence ###
  z_ijtk[] <- rbinom(n=length(c(gamma_ijtk)),size=1,prob=c(lambda_ijtk))
  ### link weight, assuming there's a link ###
  y_ijtk[] <- rnorm(n=length(c(mu_ijtk)),mean=c(mu_ijtk),sd=rep(sigma_k,each=V_net*V_net*T_net))
  
  low_tri <- lower.tri(y_ijtk[,,1,1])
  
  for( k in 1:K_net ) { #k<-1
    for( t in 1:T_net ) { #i<-1
      diag(z_ijtk[,,t,k]) <- diag(y_ijtk[,,t,k]) <- 0
    }
  }
  if(!directed){
    for( t in 1:T_net) {
      for( k in 1:K_net) { #t<-1;k<-1
        z_ijtk[,,t,k][t(low_tri)] <- z_ijtk[,,t,k][low_tri]
        y_ijtk[,,t,k][t(low_tri)] <- y_ijtk[,,t,k][low_tri]
      }
    }
  }
  
  ### Weights ###
  if(!weighted){
    y_ijtk = z_ijtk
  } else {
    y_ijtk = z_ijtk * y_ijtk
  }
  
  
  ### Covariates ###
  P_pred=2
  aux <- cbind( baseline_1(times_net,0.3*T_net,0.2*T_net,3,-1.5,-1),
                baseline_1(times_net,0.3*T_net,0.2*T_net,0.5,0,1) )
  # plot(y=aux[,1],x=times_net)
  # plot(y=aux[,2],x=times_net)
  x_ijtkp <- array( data=NA, dim=c(V_net,V_net,T_net,K_net,P_pred) )
  x_ijtkp[,,,,1] <- rnorm(n=length(c(x_ijtkp[,,,,1])),mean=0,sd=0.05)
  x_ijtkp[,,,,2] <- rnorm(n=length(c(x_ijtkp[,,,,1])),mean=0,sd=0.03)
  for( l in 1:P_pred) { #l <- 1
    for( t in 1:T_net) { #t <- 1
      x_ijtkp[,,t,,l] <- x_ijtkp[,,t,,l] + aux[t,l]
    }
  }
  
  synth_net <- list( y_ijtk = exp(y_ijtk)-1,
                     
                     x_ijtkp = x_ijtkp,
                     
                     nodes_net=nodes_net,
                     times_net=times_net,
                     layers_net=layers_net,
                     
                     directed=directed,
                     weighted=weighted,
                     
                     mu_ijtk = mu_ijtk,
                     lambda_ijtk = lambda_ijtk,
                     gamma_ijtk = gamma_ijtk )
  
  if(!is.null(rds_file)){
    if( substr(rds_file,nchar(rds_file)-3,nchar(rds_file))!=".rds" ) {rds_file<-paste(rds_file,".rds",sep="")}
    saveRDS( synth_net, file=rds_file )
  }
  return( synth_net )
}
