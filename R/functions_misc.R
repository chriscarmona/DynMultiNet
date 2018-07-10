#' @importFrom stats var
R_hat.mcmc <- function(x,m) {
  # Computes the Potential Scale Reduction Coefficient
  # Gelman et al. (2014) sec. 11.4 page 285
  n <- floor(length(x)/m)
  chain_split <- matrix(NA,nrow=n,ncol=m)
  for(i in 1:m) {
    chain_split[,i] <- x[(1+n*(i-1)):(n*i)]
  }
  B <- n * stats::var(apply(chain_split,2,mean))
  W <- mean(apply(chain_split,2,stats::var))
  var_hat <- ((n-1)/n)*W + (1/n)*B
  R_hat <- sqrt(var_hat/W)
  R_hat
}

#' @export
bounce_limit <- function(x,a,b){
  while( (x<a) || (x>b) ) {
    if(x < a) {
      x <- a + (a-x)
    }
    if(x > b) {
      x <- b - (x-b)
    }
  }
  return(x)
}

#' @import grid
#' @export
multiplot <- function(..., plotlist=NULL, cols) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    grid::viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
  
}