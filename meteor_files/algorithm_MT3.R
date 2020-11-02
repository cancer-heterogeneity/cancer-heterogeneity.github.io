program <-
  function(D_matrix, k=5) {
    ##
    ## YOUR CODE BEGINS HERE
    ##
    
    if ( !{ "medepir" %in% installed.packages( ) } ) {
      devtools::install_github(repo = "bcm-uga/medepir")
    }
    if ( !{ "Edec" %in% installed.packages( ) } ) {
      devtools::install_github(repo = "BRL-BCM/EDec")
    }
    library(medepir)
    A_matrix <- NULL
    if (!is.null(x = D_matrix) ) {
      
      D_met = medepir::feature_selection(D_matrix)
      res = medepir::Edec(D_met, nbcell =k, infloci = rownames(D_met))
      A_matrix <-  res$A
      T_matrix <- res$T
      remove(list = "res")
    }
    
    ##
    ## YOUR CODE ENDS HERE
    ##
    
    return( list(A_matrix = A_matrix,T_matrix = T_matrix) )
    
  }
