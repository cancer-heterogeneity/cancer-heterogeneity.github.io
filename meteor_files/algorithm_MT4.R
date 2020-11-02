program <-
  function(D_matrix, k=5) {
    ##
    ## YOUR CODE BEGINS HERE
    ##
    
    if ( !{ "medepir" %in% installed.packages( ) } ) {
      devtools::install_github(repo = "bcm-uga/medepir")
    }
    if ( !{ "MeDeCom" %in% installed.packages( ) } ) {
      devtools:::install_github(repo = "lutsik/MeDeCom")
    }
    library(package = "parallel")
    
    
    ## we compute the estimation of A for the data set :
    A_matrix <- NULL
    if ( !is.null(x = D_matrix) ) {
      D_met = medepir::feature_selection(D_met)
      res <- MeDeCom::runMeDeCom(data = D_met,
                                 Ks = k,
                                 lambda = 0.01,
                                 NCORES = 2
      )
      A_matrix <- MeDeCom::getProportions(res)
      lmcs<-MeDeCom::getLMCs(res, K=k, lambda=0.01)
      T_matrix <-  lmcs
      remove(list = "res")
    }
    
    ##
    ## YOUR CODE ENDS HERE
    ##
    
    return( list(A_matrix = A_matrix,T_matrix = T_matrix) )
    
  }
