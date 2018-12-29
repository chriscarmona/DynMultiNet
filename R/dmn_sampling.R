#' @title
#'    Bayesian Learning of Dynamic Multilayer Networks
#'
#' @description
#'    \code{dmn_sampling} Implements model from Carmona and Martinez-Jaramillo, 2018
#'
#' @param y_ijtk Array. Network information, should be dimension (V_net,V_net,T_net,K_net).
#' @param pred_data Data frame. Linked predictors information.
#' @param directed Boolean. Indicates if the provided network is directed, i.e. the adjacency matrix is assymetrical.
#' @param weighted Boolean. Indicates if the provided network is weighted, i.e. edges with values other that 0 and 1.
#' @param H_dim Integer. Latent space dimension.
#' @param R_dim Integer. Latent space dimension, for layer specific latent vectors.
#' @param add_eff_weight Boolean. Indicates if dynamic additive effects by node should be considered for edge weights.
#' @param add_eff_link Boolean. Indicates if dynamic additive effects by node should be considered for links.
#' @param delta Positive scalar. Hyperparameter controlling for the smoothness in the dynamic of latent coordinates. Larger=smoother.
#' @param shrink_lat_space Boolean. Indicates if the space should be shrinked probabilistically.
#' @param a_1 Positive scalar. Hyperparameter controlling for number of effective dimensions in the latent space.
#' @param a_2 Positive scalar. Hyperparameter controlling for number of effective dimensions in the latent space.
#' @param procrustes_lat Boolean. Indicates if the latent coordinates should be stabilised using the procustres transformation.
#' @param n_chains_mcmc Integer. Number of parallel MCMC chains.
#' @param n_iter_mcmc Integer. Number of iterations for the MCMC.
#' @param n_burn Integer. Number of iterations discarded as part of the MCMC warming up period at the beginning of the chain.
#' @param n_thin Integer. Number of iterations discarded for thining the chain (reducing the autocorrelation). We keep 1 of every n_thin iterations.
#' @param keep_y_ijtk_imp Boolean. Indicates wheter the chain with imputed missing links will be saved (FALSE by default)
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
#'                                  pred_data = NULL,
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

dmn_sampling <- function( y_ijtk,
                          pred_data=NULL,
                          directed=FALSE, weighted=FALSE,
                          
                          H_dim=10, R_dim=10,
                          
                          add_eff_weight=FALSE,
                          add_eff_link=FALSE,
                          
                          delta=36,
                          
                          shrink_lat_space=TRUE,
                          a_1=2, a_2=2.5,
                          
                          procrustes_lat=FALSE,
                          
                          n_chains_mcmc=1,
                          n_iter_mcmc=10000, n_burn=floor(n_iter_mcmc/4), n_thin=3,
                          
                          keep_y_ijtk_imp = FALSE,
                          
                          rds_file=NULL, log_file=NULL,
                          
                          quiet_mcmc=FALSE,
                          parallel_mcmc=FALSE ) {
  
  mcmc_clock <- Sys.time()
  
  # if(!is.null(pred_data)) {
  #   if( !all( is.element(unique(pred_data[,"layer"]),c(NA,unique(net_data$layer))) ) ) {
  #     stop('Layers in "pred_data" must be one of layers in "net_data"')
  #   }
  # }
  if( !is.null(rds_file) & substr(rds_file,nchar(rds_file)-3,nchar(rds_file))!=".rds" ) {
    rds_file<-paste(rds_file,".rds",sep="")
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
  
  ### Predictors data ###
  pred_net <- get_z_pred( pred_data,
                          node_all, time_all, layer_all,
                          quiet=FALSE )
  
  pred_all <- pred_net$pred_all
  
  pred_id_layer <- pred_net$pred_id_layer
  pred_id_edge <- pred_net$pred_id_edge
  
  z_tkp<-pred_net$z_tkp
  z_ijtkp<-pred_net$z_ijtkp
  
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
        "----- Inferential parameters -----\n",
        "H_dim = ",H_dim,"\n",
        "R_dim = ",R_dim,"\n",
        "delta = ",delta,"\n",
        "shrink_lat_space = ",shrink_lat_space,"\n",
        "a_1 = ",a_1,"\n",
        "a_2 = ",a_2,"\n",
        "procrustes_lat = ",procrustes_lat,"\n",
        
        "----- MCMC parameters -----\n",
        "n_iter_mcmc = ",n_iter_mcmc,"\n",
        "n_burn = ",n_burn,"\n",
        "n_thin = ",n_thin,"\n",
        "----- Storage and processing -----\n",
        "rds_file = ",rds_file,"\n",
        "log_file = ",log_file,"\n",
        "parallel_mcmc = ",parallel_mcmc,"\n",
        "---------------------------\n\n",
        "Process starting time:\n",as.character(mcmc_clock),"\n\n",
        "---------------------------\n\n",
        "MCMC Starting time:\n",as.character(Sys.time()),"\n\n",
        "---------------------------\n\n\n",
        "iter_i , mcmc_acum_minutes , Sys.time \n",
        file=log_file )
  }
  mcmc_clock <- Sys.time()
  
  if( !directed & !weighted ) {
    dmn_mcmc <- mcmc_d_0_w_0( y_ijtk=y_ijtk,
                              node_all=node_all, time_all=time_all, layer_all=layer_all,
                              
                              pred_all=pred_all,
                              pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                              z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                              
                              H_dim=H_dim, R_dim=R_dim,
                              
                              add_eff_link=add_eff_link,
                              
                              delta=delta,
                              
                              shrink_lat_space=shrink_lat_space,
                              a_1=a_1, a_2=a_2,
                              
                              procrustes_lat=procrustes_lat,
                              
                              n_iter_mcmc=n_iter_mcmc, n_burn=n_burn, n_thin=n_thin,
                              
                              rds_file=rds_file, log_file=log_file,
                              quiet_mcmc=quiet_mcmc,
                              parallel_mcmc=parallel_mcmc )
  } else if( directed & !weighted ) {
    dmn_mcmc <- mcmc_d_1_w_0( y_ijtk=y_ijtk,
                              node_all=node_all, time_all=time_all, layer_all=layer_all,
                              
                              pred_all=pred_all,
                              pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                              z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                              
                              H_dim=H_dim, R_dim=R_dim,
                              
                              add_eff_link=add_eff_link,
                              
                              delta=delta,
                              
                              shrink_lat_space=shrink_lat_space,
                              a_1=a_1, a_2=a_2,
                              
                              procrustes_lat=procrustes_lat,
                              
                              n_iter_mcmc=n_iter_mcmc, n_burn=n_burn, n_thin=n_thin,
                              
                              rds_file=rds_file, log_file=log_file,
                              quiet_mcmc=quiet_mcmc,
                              parallel_mcmc=parallel_mcmc )
  } else if( !directed & weighted ) {
    dmn_mcmc <- mcmc_d_0_w_1( y_ijtk=y_ijtk,
                              node_all=node_all, time_all=time_all, layer_all=layer_all,
                              
                              H_dim=H_dim, R_dim=R_dim,
                              
                              add_eff_link=add_eff_link,
                              add_eff_weight=add_eff_weight,
                              
                              delta=delta,
                              
                              shrink_lat_space=shrink_lat_space,
                              a_1=a_1, a_2=a_2,
                              
                              procrustes_lat=procrustes_lat,
                              
                              n_iter_mcmc=n_iter_mcmc, n_burn=n_burn, n_thin=n_thin,
                              
                              rds_file=rds_file, log_file=log_file,
                              quiet_mcmc=quiet_mcmc,
                              parallel_mcmc=parallel_mcmc )
  } else if( directed & weighted ) {
    dmn_mcmc <- mcmc_d_1_w_1( y_ijtk=y_ijtk,
                              node_all=node_all, time_all=time_all, layer_all=layer_all,
                              
                              H_dim=H_dim, R_dim=R_dim,
                              
                              add_eff_link=add_eff_link,
                              add_eff_weight=add_eff_weight,
                              
                              delta=delta,
                              
                              shrink_lat_space=shrink_lat_space,
                              a_1=a_1, a_2=a_2,
                              
                              procrustes_lat=procrustes_lat,
                              
                              n_iter_mcmc=n_iter_mcmc, n_burn=n_burn, n_thin=n_thin,
                              
                              rds_file=rds_file, log_file=log_file,
                              quiet_mcmc=quiet_mcmc,
                              parallel_mcmc=parallel_mcmc )
  } else {
    stop( "directed=",directed, ", weighted=",weighted," not currently supported.")
  }
  
  return( dmn_mcmc )
  
}
