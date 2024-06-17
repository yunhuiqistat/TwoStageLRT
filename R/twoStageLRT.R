#' Two Stage Likelihood Ratio Test for Detection of Gene Expression Mid-parent Heterosis
#'
#' This function apply two stage likelihood ratio test (2sLRT) to multi-family RNA-seq data with formats: three count matrices for famale, male and hybrid with genes in rows, samples in columns.
#' @import ggplot2
#' @import DESeq2
#' @param hybrid_mat hybrids expression count with genes in rows, samples in columns.
#' @param female_mat female expression count with genes in rows, samples in columns.
#' @param male_mat male expression count with genes in rows, samples in columns.
#' @param disp_mat if "Pois", use Pois LRT and only output Pois results; if NULL, use dispersion estimation  output will have both NBLRT and PoisLRT results; else, use provided disp_mat for NB test, output will have both NBLRT and PoisLRT results.
#' @param normalization logical, whether or not need normalization for the raw data matrices, if FALSE, treat the normalization factors all 1.
#' @param norm_factors if normalization = TRUE and norm_factors = NULL, use DESeq normalization factors; if normalization = TRUE and norm_factors is a list contains three vectors female, male and hybrid normalization factors, use the provided norm_factors in tests.
#' @param n_fam_thres genes with at least n_fam_thres null families will be used for dispersion estimation, otherwise, will use loess prediction to get dispersion, default is 20.
#' @param gene_groups_no the number of clusters for genes, default is 3.
#' @return if disp_mat = "Pois", return Poisson LRT for heterosis pvalue and adjusted pvalue from BH adjustment; if disp_mat = NULL or provided, return both NB and Poisson LRT for heterosis pvalue and qvalue and genewise dispersion, and the estimation of mu_gi, which is the female/male/hybrid_mat_normed.
#' @export
#' @examples
#' library(TwoStageLRT)
#' data(sim_data)
#' res <- twoStageLRT(female_mat = sim_data$female_mat, male_mat = sim_data$male_mat,
#'   hybrid_mat = sim_data$hybrid_mat, disp_mat = NULL,
#'   normalization = TRUE, norm_factors = NULL,
#'   n_fam_thres = 10, gene_groups_no = 2)
#' summary(res$pval_NB)
#' summary(res$pval_Pois)
#' summary(res$DISP)


twoStageLRT <- function(hybrid_mat, female_mat, male_mat, disp_mat = NULL,
                        normalization = FALSE, norm_factors = NULL,
                        n_fam_thres = 20, gene_groups_no = 3){
  n_family <- ncol(hybrid_mat)
  n_gene <- nrow(hybrid_mat)

  ##### Step 1: Determine normalization ####
  if(normalization & is.null(norm_factors)){ # use DESeq normalization
    cat("Use DESeq estimated normalization factors. \n")
    data <- cbind(female_mat, male_mat, hybrid_mat)
    DGE <- DESeqDataSetFromMatrix(countData = data,
                                  colData = data.frame(sample = 1:ncol(data)),
                                  design = model.matrix(~c(1:ncol(data))))
    DGE <- estimateSizeFactors(DGE)
    dat_normed <- DESeq2::counts(DGE, normalized = T) # this can serve as the MLE of mu_gi under global/Hg1, just y_gi/C_i
    female_mat_normed <- dat_normed[,1:n_family]
    male_mat_normed <- dat_normed[,(n_family+1):(2*n_family)]
    hybrid_mat_normed <- dat_normed[,(2*n_family+1):(3*n_family)]
    norm_factors <- list(female_norm_factor = DGE$sizeFactor[1:n_family],
                         male_norm_factor = DGE$sizeFactor[(n_family+1):(2*n_family)],
                         hybrid_norm_factor = DGE$sizeFactor[(2*n_family+1):(3*n_family)])

  }
  else if(normalization & (!is.null(norm_factors))){
    cat("Use provided normalization factors, no estimation. \n")
    female_mat_normed <- sweep(female_mat, 2, norm_factors$female_norm_factor, "/")
    male_mat_normed <- sweep(male_mat, 2, norm_factors$male_norm_factor, "/")
    hybrid_mat_normed <- sweep(hybrid_mat, 2, norm_factors$hybrid_norm_factor, "/")
  }
  else{ # when normalization is FALSE, all norm_factors will be set to be 1
    cat("Use all 1 normalization factors, no normalization. \n")
    norm_factors <- list(female_norm_factor = rep(1, n_family),
                         male_norm_factor = rep(1, n_family),
                         hybrid_norm_factor = rep(1, n_family))
    female_mat_normed <- female_mat
    male_mat_normed <- male_mat
    hybrid_mat_normed <- hybrid_mat
  }

  # from sample-wise normalization factors to a matrix for matrix computation, each column has the same normalization factor
  female_norm_factor_mat <- matrix(norm_factors$female_norm_factor, nrow = n_gene, ncol = n_family, byrow = TRUE)
  male_norm_factor_mat <- matrix(norm_factors$male_norm_factor, nrow = n_gene, ncol = n_family, byrow = TRUE)
  hybrid_norm_factor_mat <- matrix(norm_factors$hybrid_norm_factor, nrow = n_gene, ncol = n_family, byrow = TRUE)

  # Poisson LRT for heterosis
  pval_Pois <- matrix(mapply(Hetero_Pois_LRT, female_mat, male_mat, hybrid_mat,
                             female_norm_factor_mat, male_norm_factor_mat, hybrid_norm_factor_mat),
                      nrow = n_gene, ncol = n_family)
  rownames(pval_Pois) <- rownames(hybrid_mat)
  colnames(pval_Pois) <- colnames(hybrid_mat)
  qval_Pois <- apply(pval_Pois, 2, function(x){p.adjust(x, method = "fdr")})


  ##### Step 2: Gene-wise pre-test - Poisson null family LRT ####
  # compute mid-parent (mu_g1+mu_g2)/2 globally, so the estimation are (y_g1/C_1 + y_g2/C2)/2
  MP <- (female_mat_normed + male_mat_normed) / 2
  rownames(MP) <-  rownames(hybrid_mat)
  colnames(MP) <- colnames(hybrid_mat)

  # note: all the input should be raw count matrices and the normalization factor matrices
  pval_family <- matrix(mapply(Pois_3mean_LRT, female_mat, male_mat, hybrid_mat,
                               female_norm_factor_mat, male_norm_factor_mat, hybrid_norm_factor_mat),
                        nrow = n_gene, ncol = n_family)
  rownames(pval_family) <- rownames(hybrid_mat)
  colnames(pval_family) <-colnames(hybrid_mat)
  qval_family <- apply(pval_family, 2, function(x){p.adjust(x, method = "fdr")})


  if(is.null(disp_mat)){
    ##### Step 3: Estimate dispersion matrix ####
    null_disp_ind <- qval_family > 0.05 # default is 0.05, change to 0.1 in our real data application
    # SS1: for the binary matrix, use hclust with jaccard dist to group genes together and use the noull families for the worst gene as group to estimate tagwise dispersion.
    DISP <- rep(NA, nrow(null_disp_ind)) # each gene will only have one dispersion estimation
    names(DISP) <- rownames(null_disp_ind)
    null_disp_ind_sub <- null_disp_ind[apply(null_disp_ind, 1, function(x) sum(x == 1)) >= n_fam_thres,] # keep genes with n_null_fam > n_fam_thres for clustering
    jacc_clust <- cutree(hclust(dist(null_disp_ind_sub, "binary")), gene_groups_no)
    cat("Gene clusters size: \n")
    print(table(jacc_clust))
    for (i_clust in unique(jacc_clust)) { # for each cluster, find the worst gene, use its null families as groups for dispersion estimation
      genes_clust <- rownames(null_disp_ind_sub)[jacc_clust == i_clust]
      ind_clust <- null_disp_ind_sub[genes_clust,]
      worst_gene <- genes_clust[which.max(apply(ind_clust, 1, function(x) sum(x==0)))]
      families_clust <- colnames(null_disp_ind_sub)[null_disp_ind_sub[worst_gene,] == 1]
      cat(paste("Dispersion estimation for cluster", i_clust, "the number of families is", length(families_clust),", number of genes is ", length(genes_clust)," \n"))
      counts_clust <- matrix(NA, nrow = length(genes_clust), ncol = 3 * length(families_clust))
      for(n_hybrid in 1:length(families_clust)){
        counts_clust[,(n_hybrid - 1) * 3 + 1] <- female_mat[genes_clust, families_clust[n_hybrid]]
        counts_clust[,(n_hybrid - 1) * 3 + 2] <- male_mat[genes_clust, families_clust[n_hybrid]]
        counts_clust[,(n_hybrid - 1) * 3 + 3] <- hybrid_mat[genes_clust, families_clust[n_hybrid]]
      }
      group_clust <- rep(paste("G",1:length(families_clust), sep = ""), each = 3)
      # input raw counts and the edgeR will use RLE method before estimation of dispersion
      DISP[genes_clust] <- Disp_est_null_family(counts = counts_clust, group = group_clust) # include raw counts and groups, RLE=DESeq will be applied before estimation of dispersion parameters within this function
    }
    cat("Dispersion estimation summary:\n")
    print(summary(DISP))

    #SS2: fit a loess curve to predict the dispersions for genes whose number of all-null families do not reach the minimal requirement
    loess_x <- apply(MP, 1, median)
    x <- loess_x[!is.na(DISP)]
    y <- DISP[!is.na(DISP)]
    if(sd(x) > 0){
      new_x <- loess_x[is.na(DISP)]
      loess_fit <- loess(y~x, span = 0.5) # default is 0.5, change to 0.1 in our real data application
      pred_disp <- predict(loess_fit, newdata = new_x)
      pred_disp[pred_disp < 0] <- 0
      intra <- new_x[!is.na(pred_disp) & pred_disp > 0]
      intra_disp <- pred_disp[!is.na(pred_disp) & pred_disp > 0]
      extra <- new_x[is.na(pred_disp) | pred_disp <= 0]
      extra_disp <- ifelse(extra < min(x), y[which.min(x)], y[which.max(x)])
      pred_disp[is.na(pred_disp) | pred_disp <= 0] <- extra_disp
      DISP[is.na(DISP)] <- pred_disp
    }
    dispersion_mat <- matrix(DISP, nrow = n_gene, ncol = n_family)
    # visualization of single dataset dispersion estimation vs MP
    p_single_disp_vis <- ggplot()+
      geom_point(data=data.frame(x=x, y=y),aes(x=x, y=y), color = "black") +
      geom_point(data=data.frame(x=intra, y=intra_disp), aes(x=x, y=y), color="green")+
      geom_point(data=data.frame(x=extra, y=extra_disp), aes(x=x, y=y), color="blue")+
      geom_line(data=data.frame(x=x, y=loess_fit$fitted),aes(x=x, y=y), color = "red") +
      labs(title = "Single Dataset: Estimated(black), Fitted(red), \n Predicted(green), Extrapolation(blue)", x="MP", y="Dispersion")
  }
  else if(! "character" %in% class(disp_mat)){ # use provided dispersion parameters
    dispersion_mat <- disp_mat
  }
  else{ # use Poisson LRT for MPH
    cat("Output is from Poisson LRT for heterosis. \n")
    return(list(pval_Pois = pval_Pois, qval_Pois = qval_Pois,
                pval_NB = NULL, qval_NB = NULL, DISP = 0,
                female_mat_normed = NULL, male_mat_normed = NULL,
                hybrid_mat_normed = NULL, p_single_disp_vis = NULL))
  }

  ##### Step 4: NB estimation with known dispersion matrix ####
  cat("Output is from NB LRT for heterosis. \n")
  pval_NB <- matrix(mapply(Hetero_NB_LRT, female_mat, male_mat, hybrid_mat,
                           female_norm_factor_mat, male_norm_factor_mat, hybrid_norm_factor_mat,
                           dispersion_mat), nrow = n_gene, ncol = n_family)
  rownames(pval_NB) <- rownames(hybrid_mat)
  colnames(pval_NB) <- colnames(hybrid_mat)
  qval_NB <- apply(pval_NB, 2, function(x){p.adjust(x, method = "fdr")})
  return(list(pval_NB = pval_NB, qval_NB = qval_NB, DISP = dispersion_mat[,1],
              pval_Pois = pval_Pois, qval_Pois = qval_Pois,
              female_mat_normed = female_mat_normed, male_mat_normed = male_mat_normed,
              hybrid_mat_normed = hybrid_mat_normed, p_single_disp_vis = p_single_disp_vis))
}
