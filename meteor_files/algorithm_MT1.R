program <-
  function(D_matrix, k=5) {
    ##
    ## YOUR CODE BEGINS HERE
    ##
    
    if ( !{ "deconica" %in% installed.packages( ) } ) {
      devtools::install_github(repo = "UrszulaCzerwinska/DeconICA", build_vignettes = FALSE, dependencies = TRUE)
    }
    source(file = "http://sablab.net/scripts/LibICA.r")
    library(deconica)
    
    ## we compute the estimation of A for the data set :
    A_matrix <- NULL
    if ( !is.null(x = D_matrix) ) {
      
      Rna = D_matrix
      IC = runICA(Rna, ncomp = 10, ncores = 2, ntry = 50)
      ## take stable components
      idc = apply(IC$stab, 2, mean) > 0.8
      Genes = getGenesICA(IC, alpha=0.2)
      genes = unlist(lapply(Genes[idc],function(x){c(x$pos$genes, x$neg$genes)}))
      Rna <- Rna[genes, ]
      
      ### RNA ###
      # Run FASTICA
      rna_data = Rna[!duplicated(Rna), ]
      print(dim(rna_data))
      ICA_deconv = deconica::run_fastica(rna_data, gene.names = row.names(rna_data),
                                         overdecompose = FALSE,
                                         with.names = FALSE,
                                         n.comp = k,
                                         R = TRUE
      )
      # Take list of 30 most important genes
      weighted.list_30 <- deconica::generate_markers(df = ICA_deconv,
                                                     n = 30,
                                                     return = "gene.ranked"
      )
      # Use the most important genes to weight the components score
      ICA_5_scores_weighted = deconica::get_scores(ICA_deconv$log.counts, weighted.list_30, summary = "weighted.mean", na.rm = TRUE)
      
      # Extract the proportion from the weighted scores
      print("Before problematic DeconICA function")
      tmp_rna = deconica::stacked_proportions_plot(t(ICA_5_scores_weighted))
      print("After problematic DeconICA function")
      A_rna = matrix(tmp_rna$data$value, nrow = 5, ncol = 30)
      
      A_matrix <- A_rna
      T_matrix <- ICA_deconv$S
    }
    
    ##
    ## YOUR CODE ENDS HERE
    ##
    
    return( list(A_matrix = A_matrix,T_matrix = T_matrix) )
    
  }
