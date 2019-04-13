#' @title
#'    Bayesian Learning of Dynamic Multilayer Networks
#'
#' @description
#'    \code{dmn_sampling} Implements model from Carmona and Martinez-Jaramillo, 2018
#'
#' @param y_ijtk Array. Network information, should be dimension (V_net,V_net,T_net,K_net).
#' @param directed Boolean. Indicates if the provided network is directed, i.e. the adjacency matrix is assymetrical.
#' @param weighted Boolean. Indicates if the provided network is weighted, i.e. edges with values other that 0 and 1.
#' @param x_ijtkp Array. Edge Specific external covariates.
#' @param H_dim Integer. Latent space dimension, for global latent vectors.
#' @param R_dim Integer. Latent space dimension, for layer specific latent vectors.
#' @param add_eff_weight Boolean. Indicates if dynamic additive effects by node should be considered for edge weights.
#' @param add_eff_link Boolean. Indicates if dynamic additive effects by node should be considered for links.
#' @param class_dyn character. Specifies the dynamics for latent elements: "GP" for Gaussian Processes, "nGP" for Nested Gaussian Processes.
#' @param delta Positive scalar. Hyperparameter controlling for the smoothness in the dynamic of latent coordinates. Larger=smoother, only valid for class_dyn="GP".
#' @param n_iter_mcmc Integer. Number of iterations for the MCMC.
#' @param n_burn Integer. Number of iterations discarded as part of the MCMC warming up period at the beginning of the chain.
#' @param n_thin Integer. Number of iterations discarded for thining the chain (reducing the autocorrelation). We keep 1 of every n_thin iterations.
#' @param slim_mcmc_out Boolean. Indicates if only the main components of the MCMC will be returned (TRUE by default)
#' @param rds_file String. Indicates a file (.rds) where the output will be saved.
#' @param log_file String. Indicates a file (.txt) where the log of the process will be saved.
#' @param quiet_mcmc Boolean. Indicates if silent mode is preferes, if \code{FALSE} progress update is displayed.
#' @param parallel_mcmc Boolean. Indicates if some steps in the mcmc would be processed in parallel.
#'
#' @details
#'    The model assumes a latent variable approach
#'    
#'    \code{net_data} must be a data frame with the following columns
#'        \describe{
#'        \item{\code{source}}{Start Node}
#'        \item{\code{target}}{End node}
#'        \item{\code{weight}}{Edge weight}
#'        \item{\code{time}}{time associated with the edge}
#'        \item{\code{layer}}{Layer associated with the edge}
#'        }
#'
#' @return
#'    A list with the following components:
#' \describe{
#'     \item{\code{theta_mcmc}}{Matrix with the chain of the parameters in the model.}
#' }
#'
#'
#' @examples
#' 
#' \dontrun{
#' synth_net <- gen_synth_net( node_all = seq(1,10),
#'                             time_all = seq(1,15),
#'                             layer_all = seq(1,3),
#'                             directed = FALSE,
#'                             H_dim = 3, R_dim = 3,
#'                             delta=36,
#'                             a_1 = 1.5, a_2 = 2.5 )
#'                             
#' dmn_mcmc <- dmn_sampling( net_data = synth_net$edge_data,
#'                                  directed = FALSE,
#'                                  H_dim = 10, R_dim = 5,
#'                                  delta=36,
#'                                  a_1 = 2, a_2 = 2,
#'                                  n_iter_mcmc = 3000, n_burn = 1000, n_thin = 2 )
#' }
#' 
#' @useDynLib DynMultiNet, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' 
#' @import dplyr
#' 
#' @export
#' 

dmn_VB <- function( y_ijtk,
                    directed=FALSE, weighted=FALSE,
                    
                    x_ijtkp=NULL,
                    
                    H_dim=10, R_dim=10,
                    
                    add_eff_weight=FALSE,
                    add_eff_link=FALSE,
                    
                    delta=36,
                    lat_mean=TRUE,
                    sigma_lat_mean=5,
                    
                    vb_algorithm=c("meanfield","fullrank")[2],
                    rds_file=NULL, log_file=NULL ) {
  
  mcmc_clock <- Sys.time()
  
  if( !is.null(rds_file) ) {
    if( substr(rds_file,nchar(rds_file)-3,nchar(rds_file))!=".rds" ) {
      rds_file<-paste(rds_file,".rds",sep="")
    }
  }
  #### Start: Processing data ####
  ### Network data ###
  
  if(dim(y_ijtk)[1]!=dim(y_ijtk)[2]){
    stop("y_ijtk should be an array with the same cardinality for dimensions 1 and 2")
  }
  if(length(dim(y_ijtk))==3){
    warning("y_ijtk was declared with only 3 dimensions, a single layer will be assumed")
    y_ijtk <- array(y_ijtk,dim=c(dim(y_ijtk),1))
  }
  
  V_net <- dim(y_ijtk)[1]
  T_net <- dim(y_ijtk)[3]
  K_net <- dim(y_ijtk)[4]
  
  node_all <- dimnames(y_ijtk)[[1]]
  if(is.null(node_all)){node_all<-1:V_net; dimnames(y_ijtk)[[1]]<-dimnames(y_ijtk)[[2]]<-node_all}
  
  time_all <- dimnames(y_ijtk)[[3]]
  if(is.null(time_all)){time_all<-1:T_net; dimnames(y_ijtk)[[3]]<-time_all}
  time_all <- as.numeric(time_all)
  if(any(is.na(time_all))){stop("dimnames(y_ijtk)[[3]] should be NULL or able to transform to a numeric value")}
  layer_all <- dimnames(y_ijtk)[[4]]
  if(is.null(layer_all)){layer_all<-1:K_net; dimnames(y_ijtk)[[4]]<-layer_all}
  
  #### End: Processing data ####
  
  
  if( !is.null(log_file) ) {
    model_des <- "Dynamic"
    if(K_net==1) {
      model_des <- paste(model_des," single-layer network,",collapse="")
    } else if(K_net>1) {
      model_des <- paste(model_des," multi-layer network,",collapse="")
    }
    if(directed) {
      model_des <- paste(model_des," directed",collapse="")
    } else {
      model_des <- paste(model_des," undirected",collapse="")
    }
    if(weighted) {
      model_des <- paste(model_des," weighted edges",collapse="")
    } else {
      model_des <- paste(model_des," unweighted edges",collapse="")
    }
    
    cat("**** DynMultiNet *****\n\n",
        "----- Network topology -----\n",
        "Nodes = ",V_net,"\n",
        "Layers = ",K_net,"\n",
        "Times steps = ",T_net,"\n",
        
        "Directed = ",directed,"\n",
        "Weighted = ",weighted,"\n",
        
        "----- External covariates -----\n",
        "P = ",ifelse(is.null(x_ijtkp),0,dim(x_ijtkp)[5]),"\n",
        
        "----- Inferential parameters -----\n",
        "H_dim = ",H_dim,"\n",
        "R_dim = ",R_dim,"\n",
        
        "delta = ",delta,"\n",
        
        "add_eff_weight = ",add_eff_weight,"\n",
        "add_eff_link = ",add_eff_link,"\n",
        
        "lat_mean = ",lat_mean,"\n",
        "sigma_lat_mean = ",sigma_lat_mean,"\n",
        
        "----- variational inference algorithm -----\n",
        "vb_algorithm = ",vb_algorithm,"\n",
        
        "----- Storage and processing -----\n",
        "rds_file = ",rds_file,"\n",
        "log_file = ",log_file,"\n",
        "---------------------------\n\n",
        "Process starting time:\n",as.character(mcmc_clock),"\n\n",
        "---------------------------\n\n",
        "MCMC Starting time:\n",as.character(Sys.time()),"\n\n",
        "---------------------------\n\n\n",
        file=log_file )
  }
  mcmc_clock <- Sys.time()
  
  dmn_mcmc <- VB_stan( y_ijtk=y_ijtk,
                       directed=directed, weighted=weighted,
                       
                       node_all=node_all, time_all=time_all, layer_all=layer_all,
                       
                       x_ijtkp=x_ijtkp,
                       
                       H_dim=H_dim, R_dim=R_dim,
                       
                       add_eff_link=add_eff_link,
                       add_eff_weight=add_eff_weight,
                       
                       delta=delta,
                       lat_mean=lat_mean,
                       sigma_lat_mean=sigma_lat_mean,
                       
                       vb_algorithm=vb_algorithm,
                       rds_file=rds_file, log_file=log_file )
  
  if( !is.null(log_file) ) {
    cat("MCMC Finish time:\n",as.character(Sys.time()),"\n\n",
        "---------------------------\n\n\n",
        file=paste(log_file,".txt",sep=""),append=T )
  }
  
  return( dmn_mcmc )
  
}
