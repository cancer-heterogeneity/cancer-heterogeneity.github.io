program <-
function(D_matrix, k=NULL, simcan=NULL, cancer_type) {
  ##
  ## YOUR CODE BEGINS HERE
  ##
  
  if ( !{ "immunedeconv" %in% installed.packages( ) } ) {
    if (!requireNamespace("remotes", quietly = TRUE))
      install.packages("remotes")
    remotes::install_github("icbi-lab/immunedeconv")
  }
  
  library(immunedeconv)
  
  # Provide CIBERSORT files
  set_cibersort_binary(paste0(input,"/input/CIBERSORT.R"))
  set_cibersort_mat(paste0(input,"/input/LM22.txt"))
  
  ## we compute the estimation of A for the data set :
  A_matrix <- NULL
  if ( !is.null(x = D_matrix) ) {
    ## temporary hack : we keep only the 10 000 first lines in order to save time for the baseline
    if ( nrow(x = D_matrix) > 1e4 ) {
      D_matrix[seq_len(length.out = 1e4), ]
    }
    
    res <- immunedeconv::deconvolute(gene_expression = D_matrix,
                                     method = "cibersort",
                                     arrays = FALSE)
    res = data.frame(res, row.names = 1)
    #rownames(res) = paste0(rownames(res),"|","CIBERSORT")
    
    A_matrix = res
    
    ref = read.csv2(paste0(input,"/input/LM22.txt"), sep="\t")
    T_matrix = data.frame(ref, row.names=1)
  }
  
  ##
  ## YOUR CODE ENDS HERE
  ##
  
  return( list(A_matrix = A_matrix, T_matrix = T_matrix) )
  
}
