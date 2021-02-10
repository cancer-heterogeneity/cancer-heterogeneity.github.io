program <-
function(D_matrix,  k = 5, cancer_type = cancer_type){
  
  colnames(D_matrix) = paste0("sample_", 1:ncol(D_matrix))
  rna_data = D_matrix

  ### DECONVOLUTION PART : 
  ### ICA, weighting of 30 most important genes on the score 
  
  rna_data = rna_data[!duplicated(rna_data), ]
  print(dim(rna_data))
  # ICA_deconv = deconica::run_fastica(rna_data, gene.names = row.names(rna_data),
  #                                 overdecompose = FALSE,
  #                                 with.names = FALSE,
  #                                 n.comp = k,
  #                                 R = TRUE)
                                    
    # Take list of 30 most important genes
  ICA_deconv = fastICA::fastICA(X = rna_data, n.comp = k, maxit = 1000, tol = 1e-09)
  ICA_deconv$names = row.names(rna_data)
  weighted.list <- deconica::generate_markers(df = ICA_deconv,
                                                 n = 30,
                                                 return = "gene.ranked"
  )
  
  #    row.names(weighted.list$IC1) = weighted.list$IC1[,1]
  #   row.names(weighted.list$IC2) = weighted.list$IC2[,1]
  #   row.names(weighted.list$IC3) = weighted.list$IC3[,1]
  #   row.names(weighted.list$IC4) = weighted.list$IC4[,1]

    # Use the most important genes to weight the components score
  ICA_scores_weighted = deconica::get_scores(ICA_deconv$X, weighted.list, summary = "weighted.mean", na.rm = TRUE)
  
  
  
    # Extract the proportion from the weighted scores
    tmp_dat = t(ICA_scores_weighted)
    colnames(tmp_dat) = colnames(D_matrix)
    # tmp_rna = deconica::stacked_proportions_plot(tmp_dat)
    A_rna = abs(tmp_dat) %*% diag(1/colSums(abs(tmp_dat)))
    #A_rna = matrix(tmp_rna$data$value, nrow = k)
        T_matrix <- ICA_deconv$S
  
  return(list(A_matrix=A_rna,T_matrix=T_matrix))
  
}

# df = ICA_deconv$log.counts
# markers.list = weighted.list
# summary = "weighted.mean"
#
# deconica::get_scores
# function (df, markers.list, summary = "mean", ...)
# {
#     type <- NULL
#     if (summary == "weighted.mean") {
#         type = "metagenes"
#     }
#     else {
#         type = "gene.list"
#     }
#     switch(type, gene.list = sapply(markers.list, function(markers) {
#         if (!is.null(ncol(markers))) markers <- markers[, 1]
#         fun <- match.fun(summary)
#         common.markers <- intersect(markers, row.names(df))
#         if (length(common.markers) < 0.5 * length(markers.list$markers)) {
#             warning(paste("not enough markers for"), markers,
#             sep = " ")
#             return(NA)
#         } else {
#             df.m <- df[common.markers, ]
#             apply(df.m, 2, fun, ...)
#         }
#     }), metagenes = sapply(markers.list, function(metagene) {
#         markers <- metagene[, 1]
#         print(metagene)
#         print(paste0("length of markers",length(markers)))
#         common.markers <- intersect(markers, row.names(df))
#         print(paste0("length of common marker",length(common.markers)))
#         if (length(common.markers) < 0.5 * length(markers.list$markers)) {
#             warning(paste("not enough markers for"), markers,
#             sep = " ")
#             return(NA)
#         } else {
#             df.m <- df[common.markers, ]
#             apply(df.m, 2, stats::weighted.mean, w = metagene[,
#             2])
#         }
#     }))
# }
#
# row.names(D_matrix)[!(row.names(D_matrix) %in% row.names(ICA_deconv$log.counts)]
# which(row.names(ICA_deconv$log.counts)  %in% (row.names(D_matrix)))
# row.names(ICA_deconv$log.counts) [-(which(row.names(ICA_deconv$log.counts)  %in% (row.names(D_matrix))))]
# row.names(ICA_deconv$log.counts) [-(which(row.names(ICA_deconv$S)  %in% (row.names(D_matrix))))]
#
# D_matrix["ZSCAN16.AS1",]
# rna_data["ZSCAN16.AS1",]
# ICA_deconv$log.counts["ZSCAN16.AS1",]
# ICA_deconv$S["ZSCAN16.AS1",]
# ICA_deconv$X["ZSCAN16.AS1",]
#
# sum(row.names(ICA_deconv$S) == row.names(ICA_deconv$log.counts))
#
# sum(row.names(D_matrix) %in% row.names(ICA_deconv$log.counts))
# "RNU4-26P" %in% row.names(D_matrix)
# dim(D_matrix)
# "RNU4-26P" %in% row.names(ICA_deconv$log.counts)
# dim(ICA_deconv$log.counts)
# "RNU4-26P" %in% row.names(ICA_deconv$X)
#
# "ZYG11AP1" %in% row.names(D_matrix)
# "ZYG11AP1" %in% row.names(ICA_deconv$log.counts)
#
#
# df = ICA_deconv$log.counts
# markers.list = weighted.list
# summary = "weighted.mean"
#
# res = fastICA::fastICA(X = D_matrix, n.comp = 4, maxit = 1000, tol = 1e-09)
#
# row.names(D_matrix)[!(row.names(D_matrix) %in% row.names(res$X)]
# which(row.names(ICA_deconv$log.counts)  %in% (row.names(D_matrix)))
#
# row.names(ICA_deconv$log.counts) [-(which(row.names(res$X)  %in% (row.names(D_matrix))))]
#
# row.names(res$X)
# row.names(res$log.counts)
# row.names(ICA_deconv$log.counts) [-(which(row.names(ICA_deconv$S)  %in% (row.names(D_matrix))))]
#
#
# df = ICA_deconv
# sel.comp = paste("IC", 1:ncol(df$S),sep = "")
# > deconica::generate_markers
# function (df, n = 30, thr = Inf, sel.comp = paste("IC", 1:ncol(df$S),
# sep = ""), return = "gene.list", orient.long = TRUE)
# {
#     if (orient.long) {
#         S_or <- .orient_funct(df$S)
#     }
#     else {
#         S_or <- df$S
#     }
#     colnames(S_or) <- paste("IC", 1:ncol(S_or), sep = "")
#     row.names(S_or) <- df$names
#     S_or <- S_or[, sel.comp]
#     metagenes <- apply(S_or, 2, function(col) data.frame(GENE = row.names(S_or),
#     col))
#     switch(return, gene.list = lapply(metagenes, function(x) {
#         x <- x[order(-x[, 2]), ]
#         x <- x[which(x[, 2] < thr), ]
#         return(as.character(x[1:n, 1]))
#     }), gene.ranked = lapply(metagenes, function(x) {
#         x <- x[order(-x[, 2]), ]
#         return(x[which(x[, 2] < thr), ][1:n, ])
#     }))
# }
#
# deconica::stacked_proportions_plot
# function (dat)
# {
#     packages <- c("ggplot2", "reshape")
#     new.packages <- packages[!(packages %in% utils::installed.packages()[,
#     "Package"])]
#     if (length(new.packages))
#     utils::install.packages(new.packages, repos = "http://cran.us.r-project.org")
#     cols <- colnames(dat)
#     dat <- as.matrix(dat)
#     dat <- dat %*% diag(1/colSums(dat))
#     dat <- data.frame(dat)
#     colnames(dat) <- cols
#     dat$row <- row.names(dat)
#     dat2 <- reshape::melt(dat, id.vars = "row")
#     colnames(dat2)[1] <- "cell_type"
#     dat2[, 1] <- as.factor(dat2[, 1])
#     variable <- value <- cell_type <- NULL
#     ggplot2::ggplot(dat2, ggplot2::aes(x = variable, y = value,
#     fill = cell_type)) + ggplot2::geom_bar(stat = "identity") +
#     ggplot2::xlab("\nSample") + ggplot2::ylab("Relative proportion\n") +
#     ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
#     hjust = 1))
# }
# <bytecode: 0x7fd5e8475d58>
# <environment: namespace:deconica>
#
#
# dat = t(ICA_scores_weighted)
# colnames(dat) = colnames(D_matrix)
# cols <- colnames(dat)
# dat <- as.matrix(dat)
# dat <- dat %*% diag(1/colSums(dat))
# dat <- data.frame(dat)
# colnames(dat) <- cols
# dat$row <- row.names(dat)
# dat2 <- reshape::melt(dat, id.vars = "row")
# colnames(dat2)[1] <- "cell_type"
# dat2[, 1] <- as.factor(dat2[, 1])
# variable <- value <- cell_type <- NULL
