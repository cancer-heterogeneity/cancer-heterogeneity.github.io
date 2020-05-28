program <-
function(input) {
    ##
    ## YOUR CODE BEGINS HERE
    ##
    
    if ( !{ "NMF" %in% installed.packages( ) } ) {
        install.packages(pkgs = "NMF", repos = "https://cloud.r-project.org")
    }

    ## we compute the estimation of A for each D matrix inside the 'input' argument :
    A_all <- lapply(
        X   = input
      , FUN = function( D ) {
          ## we do a basic nmf factorization with k = 4
          res <- NMF::nmf(x = D, rank = 4, method = "lee")
          A   <- apply(
              X      = res@fit@H
            , MARGIN = 2
            , FUN    = function( x ) {
                x / sum( x )
            }
          )
          A
      }
    )

    ## if we have computed several estimations of A, we compute the mean of all the estimations :
    if ( length(x = A_all) > 1 ) {
        ## we stock all the estimated A matrices as a array :
        A_all <- array(
            data = unlist(x = A_all)
          , dim  = c( dim( A_all[[ 1 ]] ), length(x = A_all) )
        )

        ## we then compute the mean of all the estimated A matrices as the final A matrix :
        A <- apply(
            X      = A_all
          , MARGIN = c(1, 2)
          , FUN    = mean
        )
    } else {
        ## otherwise, we send the previously estimation of A :
        A <- A_all[[ 1 ]]
    }

    return( A )
    
    ##
    ## YOUR CODE ENDS HERE
    ##
}
dataType <-
"rna"
