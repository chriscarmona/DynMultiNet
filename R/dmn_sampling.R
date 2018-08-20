#' @title
#'    Bayesian Learning of Dynamic Multilayer Networks with binary data
#'
#' @description
#'    \code{dmn_sampling} Implements model from Durante and Dunson, 2018
#'
#' @param net_data Data frame.Network information.
#' @param pred_data Data frame. Linked predictors information.
#' @param directed Boolean. Indicates if the provided network is directed, i.e. the adjacency matrix is assymetrical.
#' @param weighted Boolean. Indicates if the provided network is weighted, i.e. edges with values other that 0 and 1.
#' @param H_dim Integer. Latent space dimension.
#' @param R_dim Integer. Latent space dimension, for layer specific latent vectors.
#' @param k_x Positive scalar. Hyperparameter controlling for the smoothness in the dynamic of latent coordinates. Smaller=smoother.
#' @param k_mu Positive scalar. Hyperparameter controlling for the smoothness in the dynamic of the baseline process. Smaller=smoother.
#' @param k_p Positive scalar. Hyperparameter controlling for the smoothness in the dynamic of the predictor coefficients. Smaller=smoother.
#' @param a_1 Positive scalar. Hyperparameter controlling for number of effective dimensions in the latent space.
#' @param a_2 Positive scalar. Hyperparameter controlling for number of effective dimensions in the latent space.
#' @param n_chains_mcmc Integer. Number of chains for the MCMC.
#' @param n_iter_mcmc Integer. Number of iterations for the MCMC.
#' @param n_burn Integer. Number of iterations discarded as part of the MCMC warming up period at the beginning of the chain.
#' @param n_thin Integer. Number of iterations discarded for thining the chain (reducing the autocorrelation). We keep 1 of every n_thin iterations.
#' @param out_file String. Indicates a file (.RData) where the output will be saved.
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
#' 
#' set.seed(0)
#' 
#' synth_net <- gen_synth_net( node_all = seq(1,10),
#'                             time_all = seq(1,15),
#'                             layer_all = seq(1,3),
#'                             directed = FALSE,
#'                             H_dim = 3, R_dim = 3,
#'                             k_x = 0.10, k_mu = 0.10, k_p = 0.10,
#'                             a_1 = 1.5, a_2 = 2.5 )
#' 
#' set.seed(0)
#' net_mcmc <- dmn_sampling( net_data = synth_net$edge_data,
#'                           pred_data = NULL,
#'                           directed = FALSE,
#'                           H_dim = 10, R_dim = 5,
#'                           k_x = 0.10, k_mu = 0.10, k_p = 0.10,
#'                           a_1 = 2, a_2 = 2,
#'                           n_iter_mcmc = 2000 )
#' }
#' 
#' @useDynLib DynMultiNet
#' 
#' @import dplyr
#' 
#' @export
#' 

dmn_sampling <- function( net_data,
                          pred_data=NULL,
                          directed=FALSE, weighted=FALSE,
                          H_dim=10, R_dim=10,
                          k_x=0.10, k_mu=0.10, k_p=0.10,
                          a_1=2, a_2=2.5,
                          n_chains_mcmc=1,
                          n_iter_mcmc=10000, n_burn=n_iter_mcmc/2, n_thin=3,
                          out_file=NULL, log_file=NULL,
                          quiet_mcmc=FALSE,
                          parallel_mcmc=FALSE ) {
  
  mcmc_clock <- Sys.time()
  
  if(!is.null(pred_data)) {
    if( !all( is.element(unique(pred_data[,"layer"]),c(NA,unique(net_data$layer))) ) ) {
      stop('Layers in "pred_data" must be one of layers in "net_data"')
    }
  }
  
  #### Start: Processing data ####
  ### Network data ###
  y_ijtk <- get_y_ijtk_from_edges( net_data,
                                   directed=directed,
                                   weighted=weighted,
                                   self_edges=FALSE )
  node_all <- sort(unique(unlist(net_data[,c("source","target")])))
  V_net <- length(node_all)
  time_all <- sort(unique(unlist(net_data$time)))
  T_net <- length(time_all)
  layer_all <- sort(unique(unlist(net_data$layer)))
  K_net <- length(layer_all)
  
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
        "k_x = ",k_x,"\n",
        "k_mu = ",k_mu,"\n",
        "k_p = ",k_p,"\n",
        "a_1 = ",a_1,"\n",
        "a_2 = ",a_2,"\n",
        "----- MCMC parameters -----\n",
        "n_chains_mcmc = ",n_chains_mcmc,"\n",
        "n_iter_mcmc = ",n_iter_mcmc,"\n",
        "n_burn = ",n_burn,"\n",
        "n_thin = ",n_thin,"\n",
        "----- Storage and processing -----\n",
        "out_file = ",out_file,"\n",
        "log_file = ",log_file,"\n",
        "parallel_mcmc = ",parallel_mcmc,"\n",
        "---------------------------\n\n",
        "Process starting time:\n",as.character(mcmc_clock),"\n\n",
        "---------------------------\n\n",
        "MCMC Starting time:\n",as.character(Sys.time()),"\n\n",
        "---------------------------\n\n\n",
        file=log_file )
  }
  mcmc_clock <- Sys.time()
  
  dmn_mcmc <- mcmc_stan( y_ijtk=y_ijtk,
                         node_all=node_all, time_all=time_all, layer_all=layer_all,
                         
                         pred_all=pred_all,
                         pred_id_layer=pred_id_layer, pred_id_edge=pred_id_edge,
                         z_tkp=z_tkp, z_ijtkp=z_ijtkp,
                         
                         H_dim=H_dim, R_dim=R_dim,
                         k_x=k_x, k_mu=k_mu, k_p=k_p,
                         a_1=a_1, a_2=a_2,
                         
                         directed=directed,
                         weighted=weighted,
                         
                         n_chains_mcmc=n_chains_mcmc,
                         n_iter_mcmc=n_iter_mcmc, n_burn=n_burn, n_thin=n_thin,
                         
                         out_file=out_file,
                         quiet_mcmc=quiet_mcmc )
  
  if( !is.null(log_file) ) {
    cat("MCMC Finish time:\n",as.character(Sys.time()),"\n\n",
        "---------------------------\n\n\n",
        file=log_file,append=T )
  }
  
  return( dmn_mcmc )
  
}