program <-
  function(D_matrix, k=5, cancer_type = cancer_type) {
    ##
    ## YOUR CODE BEGINS HERE
    ##
    
    if ( !{ "NMF" %in% installed.packages( ) } ) {
      install.packages(pkgs = "NMF")
    }
    
    ## we compute the estimation of A for the data set :
    A_matrix <- NULL
    if (!is.null(x = D_matrix) ) {
        if (nrow(D_matrix) < 5000){
            print("number of features < 5000")
           D = D_matrix
        } else {
            print("number of features > 5000, variance based feature selection is applied")
            D = medepir::feature_selection(D_matrix)
        }
      
      res <- NMF::nmf(x = D, rank = k, method = "snmf/r", seed = 1)
      A   <- apply(
        X      = res@fit@H
        , MARGIN = 2
        , FUN    = function( x ) {
          x / sum( x )
        }
      )
      
      
      A_matrix <- A
      T_matrix <- res@fit@W
      remove(list = "res")
    }
    
    ##
    ## YOUR CODE ENDS HERE
    ##
    
    return( list(A_matrix = A_matrix,T_matrix = T_matrix) )
    
  }
