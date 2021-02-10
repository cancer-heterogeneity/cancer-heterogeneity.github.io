program <-
  function(D_matrix, k=5, cancer_type = cancer_type) {
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
      
      if (nrow(D_matrix) < 5000){
          print("number of features < 5000")
          results = medepir::Edec(D_matrix, infloci = rownames(D_matrix))
      } else {
          print("number of features > 5000, variance based feature selection is applied")
          D_matrix = medepir::feature_selection(D_matrix)
          results = medepir::Edec(D_matrix, infloci = rownames(D_matrix))
      }
      
      A_matrix <-  results$A
      T_matrix <- results$T
      remove(list = "res")
    }
    
    ##
    ## YOUR CODE ENDS HERE
    ##
    
    return( list(A_matrix = A_matrix,T_matrix = T_matrix) )
    
  }
