#' Simulation Function for 2sLRT, Poisson LRT and Estimation Method.
#'
#' This function simulates n_families with each family includes three counts (i=1-female, 2-male, 3-hybrid) from Negative Binomial distribution for n_genes.
#' Then 2sLRT, Poisson based LRT for MPH and Estimation method as described in the paper are applied to check mid-parent heterosis gene-family combinations.
#' The hypotheses for gene g in family h is Hgh0: mu_{gh1} + mu_{gh2} = 2mu_{gh3} vs Hgh1: otherwise.
#' The effect size under alternative is specified by effect factor - how many times more mu_gh3 is of (mu_gh1 + mu_gh2) / 2, 0 = Hgh0.
#' @import pROC
#' @import ggplot2
#' @import cowplot
#' @import DescTools
#' @import MASS
#' @param n_genes number of genes.
#' @param n_families number of families.
#' @param female_norm_factor a vector with normalization factors for n_families females.
#'   if norm_factor == NULL, norm_factor will be all 1;
#'   if provided and are not all 1, DESeq normalization factors will be used in 2sLRT test.
#'   if provided with length = 2, norm_factor will be simulated from Unif(norm_factor[1], norm_factor[2]), DESeq will also be used.
#'   if provided with length = 1, all samples share the same normalization factor, DESeq will not be used.
#' @param hybrid_norm_factor a vector with normalization factors for n_families males, same format as described in female_norm_factor.
#' @param male_norm_factor a vector with normalization factors for n_families hybrids, same format as described in female_norm_factor.
#' @param MP_DISP_func: character name of function describing the relationship between mid-parent and dispersion parameters.
#' @param mu_female: a matrix containing mean expression count of female samples for all genes and all families with nrow = n_genes, ncol = n_families.
#' @param mu_male: a matrix containing mean count of male samples for n_genes and n_families.
#'   if a vector with length = 2, mean mu (not Cmu) for all genes and all families in male sample will be simulated from Gamma(mu_male[1], mu_male[2]) for shape and rate, thus the mean of gamma distribution is shape/rate.
#'   if a vector with length = 1, all male samples share the same mean count.
#' @param effect_factor: a matrix containing how many times more mu_{g3} is of (mu_{g1} + mu_{g2})/2, for null cases, set it to be 0.
#'   if a vector with length = 1, all families share the same effect size.
#' @param true_disp: "None" for Pois LRT, TRUE for NB LRT with true disp, FALSE for NB LRT with estimated disp.
#' @param true_norm: TRUE for using true C_i in tests, FALSE for using DESeq normalization factors in tests
#' @param n_fam_thres: genes with at least n_fam_thres all-null families will be used for dispersion estimation, otherwise, will use loess prediction to get dispersion, default is 20
#' @param gene_groups_no: the number of clusters for genes, default is 3.
#' @param repetition: the number of replications for simulation
#' @param filename: if NULL, no save of any plots, if specified, will save plots.
#' @return a list including pvalue matrix and adjusted pvalue matrix (computed using BH method).
#' @export
#' @examples
#' library(TwoStageLRT)
#' mu_female_dense <- "matrix(rgamma(case$n_genes*case$n_families, shape = 30, rate=runif(case$n_genes, 0.05, 0.2)),nrow=case$n_genes, ncol=case$n_families)"
#' effect_factor_dense <- "rbind(matrix(ifelse(rbinom(case$n_genes/4*case$n_families, 1, 0.05)==1, case$sig_pos, 0),nrow = case$n_genes/4, ncol = case$n_families), matrix(ifelse(rbinom(case$n_genes/4*case$n_families, 1, 0.05)==1, case$sig_neg, 0),nrow = case$n_genes/4, ncol = case$n_families), matrix(ifelse(rbinom(case$n_genes/4*case$n_families, 1, 0.9)==1, case$sig_pos, 0),nrow = case$n_genes/4, ncol = case$n_families), matrix(ifelse(rbinom(case$n_genes/4*case$n_families, 1, 0.9)==1, case$sig_neg, 0),nrow = case$n_genes/4, ncol = case$n_families))"
#' case <- list(n_genes = 100, n_families = 150,
#'   sig_pos = 2, sig_neg = -2/3,
#'   mu_female = mu_female_dense,
#'   mu_male = "case$mu_female", # for same mu_female and mu male, set to be "case$mu_female"
#'   effect_factor = effect_factor_dense,
#'   female_norm_factor = c(0.8, 1.3), male_norm_factor = c(0.8, 1.3), hybrid_norm_factor = c(0.8, 1.3),
#'   MP_DISP_func = "inverse_func", true_disp = FALSE, true_norm = FALSE,
#'   n_fam_thres = 10, gene_groups_no = 2, repetition = 3, plot = TRUE, filename = NULL)
#' print(case)
#' if(is.character(case$mu_female)){case$mu_female <- eval(parse(text = case$mu_female))}
#' if(is.character(case$mu_male)){case$mu_male <- eval(parse(text = case$mu_male))}
#' if(is.character(case$effect_factor)){case$effect_factor <- eval(parse(text = case$effect_factor))}
#' res <- Hetero_pip_2sLRT(n_genes = case$n_genes, n_families = case$n_families,
#'   mu_female = case$mu_female, mu_male = case$mu_male,
#'   effect_factor = case$effect_factor,
#'   female_norm_factor = case$female_norm_factor,
#'   male_norm_factor = case$male_norm_factor,
#'   hybrid_norm_factor = case$hybrid_norm_factor,
#'   true_norm = case$true_norm,
#'   MP_DISP_func = case$MP_DISP_func, true_disp = case$true_disp,
#'   n_fam_thres = case$n_fam_thres,
#'   gene_groups_no = case$gene_groups_no,
#'   repetition = case$repetition, plot = case$plot,
#'   filename = case$filename)
#' print(res)




Hetero_pip_2sLRT <- function(n_genes, n_families,
                            mu_female, mu_male, effect_factor,
                            female_norm_factor, male_norm_factor, hybrid_norm_factor, true_norm = c(TRUE, FALSE),
                            MP_DISP_func, true_disp = c("None", TRUE, FALSE),
                            n_fam_thres = 20, gene_groups_no = 3, repetition = 150,
                            plot = TRUE, filename = NULL){


  # mu_female: n_genes * n_families
  if(length(mu_female) == 1){
    cat("All genes and families have the same mean female count.\n")
    mu_female <- matrix(mu_female, nrow = n_genes, ncol = n_families)}
  else if(length(mu_female) == 2){
    cat("Use simulated mean count from gamma with shape and rate for female. \n")
    mu_female <- matrix(rgamma(n_genes * n_families, shape = mu_female[1], rate = mu_female[2]), nrow = n_genes, byrow = T)}
  else if(nrow(mu_female) == n_genes & ncol(mu_female) == n_families){
    cat("Use provided mean count matrix for female. \n")}
  else{cat("Error: mu_female must be a matrix or a vector of length 2. \n")}

  # mu_male: n_genes * n_families
  if(length(mu_male) == 1){
    cat("All genes and families have the same mean male count. \n")
    mu_male <- matrix(mu_male, nrow = n_genes, ncol = n_families)}
  else if(length(mu_male) == 2){
    cat("Use simulated mean count from gamma with shape and rate for male. \n")
    mu_male <- matrix(rgamma(n_genes * n_families, shape = mu_male[1], rate = mu_male[2]), nrow = n_genes, byrow = T)}
  else if(nrow(mu_male) == n_genes & ncol(mu_male) == n_families){
    cat("Use provided mean count matrix for male. \n")}
  else{cat("Error: mu_male must be a matrix or a vector of length 2. \n")}

  # effect_factor: n_genes * n_families
  if(length(effect_factor) == 1){
    effect_factor <- matrix(effect_factor, nrow = n_genes, ncol = n_families)}
  else if(nrow(effect_factor) == n_genes & ncol(effect_factor) == n_families){
    cat("Use provided effect factor matrix for hybrid. \n")}
  else{cat("Error: effect factor must have length either 1 or n_genes. \n")}

  # mu_hybrid: n_genes * n_families
  mu_hybrid <- (effect_factor + 1)*((mu_female + mu_male)/2)

  # female_norm_factors:  vector with size = n_families
  # male_norm_factors:  vector with size = n_families
  # hybrid_norm_factors:  vector with size = n_families
  norm_generator <- function(norm_factor){
    if(is.null(norm_factor)){
      cat("Use all 1 normalization factors for all samples. \n")
      norm_factor <- rep(1, n_families)}
    else if(length(norm_factor) == 2){
      cat("Use simulated norm factors from uniform. \n")
      norm_factor <- c(runif(n_families, min = norm_factor[1], max = norm_factor[2]))}
    else if(length(norm_factor) == 1){
      norm_factor <- rep(norm_factor, n_families)}
    else{cat("Error: norm_factors must have either 2 or 1 elements. \n")}
    return(norm_factor)
  }
  female_norm_factor <- norm_generator(female_norm_factor)
  male_norm_factor <- norm_generator(male_norm_factor)
  hybrid_norm_factor <- norm_generator(hybrid_norm_factor)
  if(mean(female_norm_factor==1) == 1 & mean(male_norm_factor==1) == 1 &
     mean(hybrid_norm_factor == 1) == 1){normalization <- FALSE}
  else{normalization <- TRUE}
  # normalization indicates whether we will have estimated norm factors, but whether use estimation is decided by true_norm

  # lambda_gi = C_imu_gi
  female_lambda <- sweep(mu_female, 2, female_norm_factor, `*`)
  male_lambda <- sweep(mu_male, 2, male_norm_factor, `*`)
  hybrid_lambda <- sweep(mu_hybrid, 2, hybrid_norm_factor, `*`)

  # dispersion: n_genes * n_families
  mp_disp_func <- get(MP_DISP_func)
  mp <- (mu_female + mu_male)/2
  # gene wise dispersion
  dispersion <- matrix(sapply(apply(mp, 1, median), mp_disp_func),
                       nrow = n_genes, ncol = n_families)
  dispersion[dispersion < 0] <- 0.01
  print(summary(c(dispersion)))

  # simulate NB data
  pval_NB <- list()
  qval_NB <- list()
  DISP_est <- matrix(NA, nrow = repetition, ncol = n_genes)
  abs_MPH <- list()
  pval_Pois <- list()
  qval_Pois <- list()

  cell_generator <- function(lambda, dispersion){
    return(rnegbin(1, mu = lambda, theta = 1/dispersion))
  }

  for(i in 1:repetition){
    female_mat <- matrix(mapply(cell_generator, female_lambda, dispersion), nrow = n_genes)
    male_mat <- matrix(mapply(cell_generator, male_lambda, dispersion), nrow = n_genes)
    hybrid_mat <- matrix(mapply(cell_generator, hybrid_lambda, dispersion), nrow = n_genes)
    rownames(hybrid_mat) <- rownames(female_mat) <- rownames(male_mat) <- paste("Gene", 1:n_genes, sep = "")
    colnames(hybrid_mat) <- colnames(female_mat) <- colnames(male_mat) <- paste("Family", 1:n_families, sep = "")

    # determine use what type of dispersions in test
    if(true_disp == "None"){disp_mat <- "Pois"} # only use Poisson LRT for heterosis
    else if(true_disp){disp_mat <- dispersion} # use NB LRT with true disp and Pois LRT for heterosis
    else{disp_mat <- NULL} # use NB LRT with estimated disp and Pois LRT for heterosis

    # determine use what type of normalization factors in test
    if(true_norm){norm_factors <- list(female_norm_factor=female_norm_factor, # use true normalization factors in test if normalization = TRUE
                                       male_norm_factor=male_norm_factor,
                                       hybrid_norm_factor=hybrid_norm_factor)}
    else{norm_factors <- NULL} # use DESeq estimated norm factors if normalization = TRUE
    res <- twoStageLRT(female_mat = female_mat, male_mat = male_mat,
                       hybrid_mat = hybrid_mat, disp_mat = disp_mat,
                       normalization = normalization, norm_factors = norm_factors,
                       n_fam_thres = n_fam_thres, gene_groups_no = gene_groups_no)

    pval_NB[[i]] <- res$pval_NB
    qval_NB[[i]] <- res$qval_NB
    DISP_est[i,] <- res$DISP
    pval_Pois[[i]] <- res$pval_Pois
    qval_Pois[[i]] <- res$qval_Pois
    abs_MPH[[i]] <- abs((res$hybrid_mat_normed - (res$female_mat_normed + res$male_mat_normed)/2)/((res$female_mat_normed + res$male_mat_normed)/2))

    ## sanity check especially for NB LRT when numeric algorithm fails
    cat(paste("Proportion of Failure NB LRT is", mean(is.na(pval_NB[[i]])), "\n"))
    if(i==1){p_single_disp_vis <- res$p_single_disp_vis}
  }

  if(plot){
    # dispersion visualization, true, estimated, true with mp, estimated with MP (on mu not y)
    p_disp <- ggplot(data = data.frame(disp = c(dispersion)), aes(x = disp))+
      geom_histogram(aes(y = after_stat(count / sum(count))))+
      labs(title = "Simulated true dispersions")
    p_disp_used <- ggplot(data = data.frame(disp = apply(DISP_est, 2, mean)), aes(x = disp))+
      geom_histogram(aes(y = after_stat(count / sum(count))))+
      labs(title = "Estimated genewise dispersions \n averaged over repetitions")
    p_disp_scatter <- ggplot(data = data.frame(mp = apply(mp, 1, median), disp = dispersion[,1]))+
      geom_point(aes(x=mp, y=disp))+
      labs(title = "True MP vs Disp")
    disp_disp_data = data.frame(Dispersion = rep(sort(dispersion[,1]),1),
                                Estimation = apply(DISP_est, 2, mean)[order(dispersion[,1])] # mean across repetitions
                                # Estimation = c(apply(DISP_est, 2, mean)[order(dispersion[,1])], # mean across repetitions
                                #                apply(DISP_est, 2, median)[order(dispersion[,1])], # median across repetitions
                                #                DISP_est[1,order(dispersion[,1])]), # for the first repetition, single dataset
                                # Type = rep(c("Mean across datasets", "Median across datasets", "Single dataset"), each = nrow(dispersion))
    )
    cat("The summary over all genes of SD(dispersion estimation across repetitions for each gene) \n")
    print(summary(apply(DISP_est, 2, sd)))
    p_disp_disp <- ggplot(data = disp_disp_data)+
      geom_point(aes(x = Dispersion, y = Estimation), size=1, alpha=0.6, color = "darkgray")+
      geom_abline(aes(slope = 1, intercept = 0, color = "x=y"), linetype = "dashed",size= 1.5)+
      theme_bw()+
      theme(legend.position = "bottom")+
      labs(title = "", color = "")
    # qq plot of pval for null case (effect_factor == 0)
    if(sum(effect_factor == 0 ) > 0){
      Empirical_NB <- quantile(unlist(lapply(pval_NB, function(x){x[effect_factor == 0]})),
                               seq(0.01,1,0.01), na.rm = TRUE)
      Empirical_Pois <- quantile(unlist(lapply(pval_Pois, function(x){x[effect_factor == 0]})),
                                 seq(0.01,1,0.01), na.rm = TRUE)
      Unif <- qunif(seq(0.01,1,0.01))
      pval_ccc_nb <- round(CCC(Unif, Empirical_NB)$rho.c["est"],4)
      pval_ccc_pois <- round(CCC(Unif, Empirical_Pois)$rho.c["est"],4)
      qq0 <- ggplot(data = data.frame(Unif = rep(Unif,2), Empirical = c(Empirical_NB, Empirical_Pois),
                                      Method = rep(c("2sLRT","Poisson"), each = length(Unif))))+
        geom_point(aes(x = Unif, y = Empirical, color=Method))+
        geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
        labs(x = "Quantiles of unif(0,1)", y = "Quantiles of pvalues", title = "", color = "")+
        xlim(0, 1)+
        ylim(0, 1)+
        theme_bw()+
        theme(legend.position = "bottom")
      cat(paste("\n CCC between pvalue and Uniform quantiles, NB=", pval_ccc_nb, ", Pois=", pval_ccc_pois, sep =""))
    }
    else{
      qq0 <- NULL
      pval_ccc_nb <- NA
    }


    if(mean(effect_factor == 0) != 1){
      # average FDR and TPR curve of qvalue with nominal alpha in (0,0.2)
      alpha_func <- function(qvalues){
        alphas <- round(seq(0.005,0.20,0.001),3)
        TPR_avg <- c()
        FPR_avg <- c()
        FDR_avg <- c()
        for (a in alphas) {
          TP_all <- lapply(qvalues, function(x){sum(x[effect_factor!=0]<=a,na.rm=TRUE)})
          FP_all <- lapply(qvalues, function(x){sum(x[effect_factor==0]<=a,na.rm=TRUE)})
          TN_all <- lapply(qvalues, function(x){sum(x[effect_factor==0]>a,na.rm=TRUE)})
          FN_all <- lapply(qvalues, function(x){sum(x[effect_factor!=0]>a,na.rm=TRUE)})
          TPR_all <- unlist(TP_all) / (unlist(TP_all) + unlist(FN_all))
          TPR_all[is.infinite(TPR_all)] <- 0
          TPR_avg <- c(TPR_avg, mean(TPR_all))
          FPR_all <- unlist(FP_all) / (unlist(FP_all) + unlist(TN_all))
          FPR_all[is.infinite(FPR_all)] <- 0
          FPR_avg <- c(FPR_avg, mean(FPR_all))
          FDR_all <- unlist(FP_all) / (unlist(TP_all) + unlist(FP_all))
          FDR_all[is.infinite(FDR_all)|is.nan(FDR_all)|is.na(FDR_all)] <- 0
          FDR_avg <- c(FDR_avg, mean(FDR_all))
        }
        return(list(alphas = alphas, FDR_avg = FDR_avg, TPR_avg = TPR_avg))
      }

      alpha_res_NB <- alpha_func(qval_NB)
      alphas <- alpha_res_NB$alphas
      FDR_avg_NB <- alpha_res_NB$FDR_avg
      TPR_avg_NB <- alpha_res_NB$TPR_avg
      alpha_res_Pois <- alpha_func(qval_Pois)
      FDR_avg_Pois <- alpha_res_Pois$FDR_avg
      TPR_avg_Pois <- alpha_res_Pois$TPR_avg
      cat(paste("\n FDR at alpha=0.05 is, NB=", round(FDR_avg_NB[alphas == 0.05],4),", Pois=",  round(FDR_avg_Pois[alphas == 0.05],4), sep =""))
      p_FDR_avg <- ggplot(data = data.frame(alpha = rep(alphas,2), FDR = c(FDR_avg_NB, FDR_avg_Pois),
                                            Method = rep(c("2sLRT","Poisson"), each = length(alphas))))+
        geom_line(aes(x = alpha, y = FDR, color = Method), linewidth = 1)+
        geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
        xlim(0, 0.2)+
        labs(title = "", y="Average FDR", x="Nominal FDR", color = "")+
        theme_bw()+
        theme(legend.position = "bottom")
      cat(paste("\n TPR at alpha=0.05 is, NB=", round(TPR_avg_NB[alphas == 0.05],4),", Pois=",  round(TPR_avg_Pois[alphas == 0.05],4), sep =""))
      p_TPR_avg <- ggplot(data = data.frame(alpha = rep(alphas,2), TPR = c(TPR_avg_NB, TPR_avg_Pois),
                                            Method = rep(c("2sLRT","Poisson"), each = length(alphas))))+
        geom_line(aes(x = alpha, y = TPR, color = Method), linewidth = 1)+
        # geom_vline(xintercept = 0.05, linetype = "dashed")+
        # geom_hline(yintercept = TPR_avg_NB[alphas == 0.05], linetype = "dashed")+
        xlim(0, 0.2)+
        ylim(0, 1)+
        labs(title = "", y="Average TPR", x="Nominal FDR", color = "")+
        theme_bw()+
        theme(legend.position = "bottom")

      at_actualFDR <- function(qvalues, actual_FDR_level = 0.05){
        ## get a list of dataframes including cutoff = sorted qvalues, FDR, FPR and TPR, one df is for one repe
        emp_values_single <- function(qvalue){
          emp_df <- data.frame(cutoff = sort(qvalue), FDR = NA, TPR = NA, FPR = NA)
          emp_df <- t(apply(emp_df, 1, function(x){
            a <- x["cutoff"]
            TP <- sum(qvalue[effect_factor!=0]<=a,na.rm=TRUE)
            FP <- sum(qvalue[effect_factor==0]<=a,na.rm=TRUE)
            TN <- sum(qvalue[effect_factor==0]>a,na.rm=TRUE)
            FN <- sum(qvalue[effect_factor!=0]>a,na.rm=TRUE)
            x["FDR"] <- ifelse(TP+FP==0, 0, FP/(TP+FP))
            x["FPR"] <- ifelse(TN+FP==0, 0, FP/(TN+FP))
            x["TPR"] <- ifelse(TP+FN==0, 0, TP/(TP+FN))
            return(x)
          }))
          return(as.data.frame(emp_df))
        }
        ## find the average of TPR, for each repetition where TPR is the max TPR whose FDR <= 0.05
        return(mean(unlist(lapply(lapply(qvalues, emp_values_single),
                                  function(x){ifelse(sum(x$FDR<=actual_FDR_level) > 0,
                                                     max(x$TPR[x$FDR <= actual_FDR_level]), NA)})),
                    na.rm=TRUE))
      }

      # cat(paste("\n Average TPR_max with empirical FDR <= 0.05 is, NB=",
      #           round(at_actualFDR(qval_NB),4),
      #           ", Pois=",  round(at_actualFDR(qval_Pois),4), "\n", sep =""))

      # average and SE AUC and partial AUC only for FPR in 0,0.2 and mean ROC curve
      mean_ROC_func <- function(qvalues, direction = c(">", "<")){
        AUC <- unlist(lapply(qvalues, function(x){roc(response = ifelse(c(effect_factor) != 0, 1, 0), predictor = c(x), direction=direction, plot = FALSE)$auc}))
        AUC_avg <- round(mean(AUC, na.rm = TRUE),4)
        AUC_SE <- round(sd(AUC, na.rm = TRUE), 4)

        partial_AUC <- unlist(lapply(qvalues, function(x){auc(roc(response = ifelse(c(effect_factor) != 0, 1, 0), predictor = c(x), direction=direction, plot = FALSE), partial.auc = c(0.8, 1), partial.auc.focus = "specificity", partial.auc.correct=TRUE)}))
        partial_AUC_avg <- round(mean(partial_AUC, na.rm = TRUE),4)
        partial_AUC_SE <- round(sd(partial_AUC, na.rm = TRUE), 4)

        roc_res <- lapply(qvalues, function(x){do.call(cbind, roc(response = ifelse(c(effect_factor) != 0, 1, 0), predictor = c(x), direction=direction, plot = FALSE)[c("sensitivities", "specificities")])})
        ## set actual FPR cutoffs (x-axis of mean ROC curve)
        fpr_actual <- seq(0.0001, 0.2, 0.001)
        ## for each actual FPR cutoff, each dataset, find the max TPR with FPR <= actual FPR
        actual_res <- lapply(roc_res, function(x){
          tpr_roc <- x[,"sensitivities"]
          fpr_roc <- 1 - x[, "specificities"]
          return(sapply(fpr_actual, function(x){max(tpr_roc[fpr_roc <= x])}))
        })
        ## for each actual FPR, compute average TPR across the datasets as the y-axis of mean ROC
        tpr_actual <- apply(do.call(rbind, actual_res), 2, mean)
        return(list(AUC_avg = AUC_avg, AUC_SE = AUC_SE,
                    partial_AUC_avg = partial_AUC_avg, partial_AUC_SE = partial_AUC_SE,
                    fpr_actual = fpr_actual, tpr_actual = tpr_actual))
      }

      NB_roc <- mean_ROC_func(qval_NB, direction = ">")
      Pois_roc <- mean_ROC_func(qval_Pois, direction = ">")
      Naive_roc <- mean_ROC_func(abs_MPH, direction = "<")

      cat(paste("\n NB: average AUC =", NB_roc$AUC_avg, " SE =", NB_roc$AUC_SE, "\n"))
      cat(paste("Pois: average AUC =", Pois_roc$AUC_avg, " SE =", Pois_roc$AUC_SE, "\n"))
      cat(paste("Naive: average AUC =", Naive_roc$AUC_avg, " SE =", Naive_roc$AUC_SE, "\n"))

      cat(paste("\n NB: average partial AUC =", NB_roc$partial_AUC_avg, " SE =", NB_roc$partial_AUC_SE, "\n"))
      cat(paste("Pois: average partial AUC =", Pois_roc$partial_AUC_avg, " SE =", Pois_roc$partial_AUC_SE, "\n"))
      cat(paste("Naive: average partial AUC =", Naive_roc$partial_AUC_avg, " SE =", Naive_roc$partial_AUC_SE, "\n"))

      p_roc_avg <- ggplot(data = data.frame(TPR = c(NB_roc$tpr_actual, Pois_roc$tpr_actual, Naive_roc$tpr_actual),
                                            FPR = rep(NB_roc$fpr_actual,3),
                                            Method = factor(rep(c("2sLRT","Poisson","Estimation"),
                                                                each = length(NB_roc$fpr_actual)),
                                                            levels=c("2sLRT","Poisson","Estimation"))))+
        geom_line(aes(x = FPR, y = TPR, color = Method, linetype = Method), linewidth = 1)+
        xlim(0, 0.2)+
        ylim(0, 1)+
        labs(title = "", color="", linetype = "")+
        theme_bw()+
        theme(legend.position = "bottom")

    }
    else{
      p_roc_avg <- NULL
      p_FDR_avg <- NULL
      p_TPR_avg <- NULL
      AUC_avg <- NA
    }
    if(is.null(qq0)){finalp <- plot_grid(p_disp, p_disp_used, p_disp_scatter, p_single_disp_vis,
                                         p_roc_avg, p_FDR_avg, p_TPR_avg, nrow = 2, ncol = 4)}
    else if(is.null(p_roc_avg)){finalp <- plot_grid(p_disp, p_disp_used, p_disp_scatter, p_single_disp_vis, qq0,
                                                    nrow = 2)}
    else{
      finalp <- plot_grid(p_disp, p_disp_used, p_disp_scatter, p_single_disp_vis,
                          qq0, p_roc_avg, p_FDR_avg, p_TPR_avg, nrow = 2, ncol = 4)
      paperp <- plot_grid(qq0, p_FDR_avg, p_roc_avg,  p_disp_disp, nrow = 1)
    }
  }

  ########### Use only on slurm #############
  if(!is.null(filename)){
    ggsave(paste(filename, "_0.pdf", sep = ""), finalp, width = 16, height = 8, units = "in")
    ggsave(paste(filename, "_4.pdf", sep = ""), paperp, width = 16, height = 4, units = "in")
  }
  ########################

  return(finalp)
}
