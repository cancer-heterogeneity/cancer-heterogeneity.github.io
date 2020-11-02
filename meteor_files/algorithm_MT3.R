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
    
    D_met = medepir::feature_selection(D_matrix)
    res = medepir::Edec(D_met, infloci = rownames(D_met))
    A_matrix <-  res$A
    T_matrix <- res$T
    remove(list = "res")
  }

##
## YOUR CODE ENDS HERE
##

return( list(A_matrix = A_matrix,T_matrix = T_matrix) )

}
