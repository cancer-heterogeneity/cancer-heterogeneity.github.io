program <-
  function(D_matrix, k=5, cancer_type) {
    ##
    ## YOUR CODE BEGINS HERE
    ##
    
    for ( package in c('NMF','csSAM','limSolve','corpcor') ) {
        if ( !{ package %in% installed.packages( ) } ) {
            install.packages(pkgs = package, repos = "https://cloud.r-project.org")
        }
    }
    
    remove(list = "package")
    
    if ( !{ 'BiocInstaller' %in% installed.packages( ) } ) {
        install.packages("BiocInstaller", repos="http://bioconductor.org/packages/3.7/bioc/")
          }
    
    # Special instructions for CellMix and DSA
     if ( !{ 'CellMix' %in% installed.packages( ) } ) {
         BiocInstaller::biocLite('CellMix', siteRepos = 'http://web.cbio.uct.ac.za/~renaud/CRAN', type='both')
     }
     
     #if ( !{ 'DSA' %in% installed.packages( ) } ) {
     #system('wget https://github.com/zhandong/DSA/raw/master/Package/version_1.0/DSA_1.0.tar.gz')
     #system("R CMD INSTALL DSA_1.0.tar.gz")
     # }
    
    
    #We generate the marker file
     mt = readRDS(paste0(input,"/marker_table.rds"))  
    mt = mt[mt$tcga_ref == cancer_type, ]
    
    #remove duplicated markers
    n_occur <- data.frame(table(mt$gene))
    non_unique = n_occur[n_occur$Freq > 1,]$Var1
    idnu = which(mt$gene %in% non_unique)
    mt = mt[-idnu, ]
    rownames(mt) = mt$gene
    
    print(paste0(nrow(mt), " genes are used as reference markers"))
    
    ## we compute the estimation of A for the data set :
    A_matrix <- NULL
    if ( !is.null(x = D_matrix) ) {
      
      T = D_matrix
      if (sum(which(rowSums(D_matrix) == 0)) > 0) {
          print("remove rows = 0")
          T = T[-which(rowSums(D_matrix) == 0),]
      }
      keep = intersect(rownames(mt),rownames(T))
      mt = mt[keep,]
      T = T[keep,]
      
      ML = CellMix::MarkerList()
      ML@.Data <- tapply(as.character(mt$gene),as.character(mt$CT),list)
      
      RESULTS <- CellMix::ged(as.matrix(T), ML, method = "ssFrobenius", sscale = FALSE, maxIter=500, log = FALSE)

      
      A_matrix <- RESULTS@fit@H
      T_matrix <- RESULTS@fit@W
    }
    
    ##
    ## YOUR CODE ENDS HERE
    ##
    
    return( list(A_matrix = A_matrix,T_matrix = T_matrix) )
    
  }
