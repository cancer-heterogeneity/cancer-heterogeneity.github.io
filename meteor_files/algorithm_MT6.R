program <-
function(D_matrix, k=NULL) {
  ##
  ## YOUR CODE BEGINS HERE
  ##
  
  if ( !{ "EpiDISH" %in% installed.packages( ) } ) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("EpiDISH")
  }
  
  
  ## we compute the estimation of A for the data set :
  A_matrix <- NULL
  if ( !is.null(x = D_matrix) ) {
    ## temporary hack : we keep only the 10 000 first lines in order to save time for the baseline
    if ( nrow(x = D_matrix) > 1e4 ) {
      D_matrix[seq_len(length.out = 1e4), ]
    }
    library(EpiDISH)
    data(centEpiFibIC.m)
    data(centBloodSub.m)
    res <- hepidish(beta.m = D_matrix, ref1.m = centEpiFibIC.m, ref2.m = centBloodSub.m, h.CT.idx = 3, method = 'RPC')
    A_matrix <- t(res)
    T_matrix = NULL
    remove(list = "res")
  }
  
  ##
  ## YOUR CODE ENDS HERE
  ##
  
  return( list(A_matrix = A_matrix,T_matrix = T_matrix) )
  
}
