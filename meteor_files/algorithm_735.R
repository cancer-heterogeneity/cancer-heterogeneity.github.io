program <-
  function(D_matrix, k=5, cancer_type = cancer_type) {
    ##
    ## YOUR CODE BEGINS HERE
    ##
    
    if ( !{ "deconica" %in% installed.packages( ) } ) {
      devtools::install_github(repo = "UrszulaCzerwinska/DeconICA", build_vignettes = FALSE, dependencies = TRUE)
    }
    
    source(paste0(input,"/input/LibICA.R"))
    if (is.null(colnames(D_matrix))) {colnames(D_matrix) = paste0("sample_", 1:ncol(D_matrix))}
    
    ## we compute the estimation of A for the data set :
    A_matrix <- NULL
    if ( !is.null(x = D_matrix) ) {
      
      print("feature selection...")
      Rna = D_matrix
      IC = runICA(Rna, ncomp = 10, ncores = 2, ntry = 50)
      ## take stable components
      idc = apply(IC$stab, 2, mean) > 0.8
      Genes = getGenesICA(IC, alpha=0.2)
      genes = unlist(lapply(Genes[idc],function(x){c(x$pos$genes, x$neg$genes)}))
      Rna <- Rna[genes, ]
     print(paste0(nrow(Rna), " features selected..."))
      
      ### RNA ###
      # Run FASTICA
      
      rna_data = Rna[!duplicated(Rna), ]
      print(dim(rna_data))
    
      ICA_deconv = fastICA::fastICA(X = rna_data, n.comp = k, maxit = 1000, tol = 1e-09)
      ICA_deconv$names = row.names(rna_data)
      weighted.list <- deconica::generate_markers(df = ICA_deconv,
      n = 30,
      return = "gene.ranked"
      )
      
      # Use the most important genes to weight the components score
      ICA_scores_weighted = deconica::get_scores(ICA_deconv$X, weighted.list, summary = "weighted.mean", na.rm = TRUE)
      
      
      
      # Extract the proportion from the weighted scores
      tmp_dat = t(ICA_scores_weighted)
      colnames(tmp_dat) = colnames(D_matrix)
      # tmp_rna = deconica::stacked_proportions_plot(tmp_dat)
     A_rna = abs(tmp_dat) %*% diag(1/colSums(abs(tmp_dat)))
     #A_rna = matrix(tmp_rna)
      
      A_matrix <- A_rna
      T_matrix <- ICA_deconv$S
    }
    
    ##
    ## YOUR CODE ENDS HERE
    ##
    
    return( list(A_matrix = A_matrix,T_matrix = T_matrix) )
    
  }

