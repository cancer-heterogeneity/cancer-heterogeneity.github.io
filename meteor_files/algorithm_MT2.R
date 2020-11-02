program <-
  function(D_matrix, k=NULL) {
    ##
    ## YOUR CODE BEGINS HERE
    ##
    
    if ( !{ "NMF" %in% installed.packages( ) } ) {
      install.packages(pkgs = "NMF")
    }
    source(file = "http://sablab.net/scripts/LibICA.r")
    
    
    ## we compute the estimation of A for the data set :
    A_matrix <- NULL
    if (!is.null(x = D_matrix) ) {
      D = D_matrix
      IC = runICA(D, ncomp=10, ncores=2, ntry=50)
      idc = apply(IC$stab, 2, mean) > 0.8
      Genes = getGenesICA(IC,alpha=0.2)
      genes = unlist(lapply(Genes[idc],function(x){c(x$pos$genes, x$neg$genes)}))
      D <- D[genes, ]
      
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
