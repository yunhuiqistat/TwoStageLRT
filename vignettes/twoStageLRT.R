library(ggplot2)
library(pROC)
library(limma)
library(MASS)
library(cowplot)
library(nleqslv)
library(rootSolve, lib.loc = "/rafalab/yunhui/Rpackages/")
library(pracma)
library(edgeR)
library(DESeq2)

args <- commandArgs(trailingOnly = TRUE)

pattern_arg <- args[1]
N0_arg <- as.numeric(args[2])
K_arg <- as.numeric(args[3])
ngenes_arg <- as.numeric(args[4])
nfamilies_arg <- as.numeric(args[5])
disp_arg <- args[6]
output_filename_arg <- args[7]

# -------------- Functions definition -------------
# This is the function CCC in DescTools
CCC <- function (x, y, ci = "z-transform", conf.level = 0.95, na.rm = FALSE)
{
    dat <- data.frame(x, y)
    if (na.rm)
        dat <- na.omit(dat)
    N. <- 1 - ((1 - conf.level)/2)
    zv <- qnorm(N., mean = 0, sd = 1)
    lower <- "lwr.ci"
    upper <- "upr.ci"
    k <- length(dat$y)
    yb <- mean(dat$y)
    sy2 <- var(dat$y) * (k - 1)/k
    sd1 <- sd(dat$y)
    xb <- mean(dat$x)
    sx2 <- var(dat$x) * (k - 1)/k
    sd2 <- sd(dat$x)
    r <- cor(dat$x, dat$y)
    sl <- r * sd1/sd2
    sxy <- r * sqrt(sx2 * sy2)
    p <- 2 * sxy/(sx2 + sy2 + (yb - xb)^2)
    delta <- (dat$x - dat$y)
    rmean <- apply(dat, MARGIN = 1, FUN = mean)
    blalt <- data.frame(mean = rmean, delta)
    v <- sd1/sd2
    u <- (yb - xb)/((sx2 * sy2)^0.25)
    C.b <- p/r
    sep = sqrt(((1 - ((r)^2)) * (p)^2 * (1 - ((p)^2))/(r)^2 +
                    (2 * (p)^3 * (1 - p) * (u)^2/r) - 0.5 * (p)^4 * (u)^4/(r)^2)/(k -
                                                                                      2))
    ll = p - zv * sep
    ul = p + zv * sep
    t <- log((1 + p)/(1 - p))/2
    set = sep/(1 - ((p)^2))
    llt = t - zv * set
    ult = t + zv * set
    llt = (exp(2 * llt) - 1)/(exp(2 * llt) + 1)
    ult = (exp(2 * ult) - 1)/(exp(2 * ult) + 1)
    if (ci == "asymptotic") {
        rho.c <- as.data.frame(cbind(p, ll, ul))
        names(rho.c) <- c("est", lower, upper)
        rval <- list(rho.c = rho.c, s.shift = v, l.shift = u,
                     C.b = C.b, blalt = blalt)
    }
    else if (ci == "z-transform") {
        rho.c <- as.data.frame(cbind(p, llt, ult))
        names(rho.c) <- c("est", lower, upper)
        rval <- list(rho.c = rho.c, s.shift = v, l.shift = u,
                     C.b = C.b, blalt = blalt)
    }
    return(rval)
}

Hetero_Pois_LRT <- function(y_g1, y_g2, y_g3, C_1 = 1, C_2 = 2, C_3 = 1){
  # This function conducts Poisson LRT test for heterosis in a single family and for a single gene with single replication per genotype
  # Hg0: mu_g1 + mu_g2 = 2 * mu_g3
  # @param y_g1: raw count for female
  # @param y_g2: raw count for male
  # @param y_g3: raw count for hybrid
  # @param C_1: normalization factor for female sample
  # @param C_2: normalization factor for male sample
  # @param C_3: normalization factor for hybrid sample
  
  # 1. Get MLE under Hg0
  if(C_1 > C_2){ # we solve mu_g2 first and then mu_g1
    a <- C_1 - C_2
    b <- 2*C_2 + C_3
    A <- a*b
    B <- b*(y_g1 + y_g2) - 2*a*(y_g2 + y_g3)
    C <- -2*y_g1*y_g2 - 2*y_g2^2 - 2*y_g2*y_g3
    # take the positive root
    mu_g2_hat0 <- ifelse((-B+sqrt(B^2 - 4*A*C)) / (2*A) > 0,
                         (-B+sqrt(B^2 - 4*A*C)) / (2*A) ,(-B-sqrt(B^2 - 4*A*C)) / (2*A))
    mu_g1_hat0 <- y_g1*mu_g2_hat0 / (y_g2 + a*mu_g2_hat0)
  }else if(C_1 < C_2){ # we solve mu_g1 first and then mu_g2
    a <- C_1 - C_2
    b <- 2*C_1 + C_3
    A <- a*b
    B <- -b*(y_g1 + y_g2) - 2*a*(y_g1 + y_g3)
    C <- 2*y_g1*y_g2 + 2*y_g1^2 + 2*y_g1*y_g3
    # take the positive root
    mu_g1_hat0 <- ifelse((-B-sqrt(B^2 - 4*A*C)) / (2*A) > 0,
                         (-B-sqrt(B^2 - 4*A*C)) / (2*A), (-B+sqrt(B^2 - 4*A*C)) / (2*A))
    mu_g2_hat0 <- y_g2*mu_g1_hat0 / (y_g1 - a*mu_g1_hat0)
  }else{ # when C_1 = C_2
    # if y_g1 + y_g2 = 0, this means y_g1=y_g2=0, means numerators are all 0, thus y_g1+y_g2 in denominator replaced 1 does not affect the MLE of mu_g12 
    mu_g1_hat0 <- 2*(y_g1^2 + y_g1*y_g2 + y_g1*y_g3) / ((2*C_1 + C_3)*ifelse(y_g1 + y_g2 == 0, 1, y_g1 + y_g2))
    mu_g2_hat0 <- 2*(y_g2^2 + y_g1*y_g2 + y_g2*y_g3) / ((2*C_2 + C_3)*ifelse(y_g1 + y_g2 == 0, 1, y_g1 + y_g2))
  }
  
  mu_g3_hat0 <- (mu_g1_hat0 + mu_g2_hat0) / 2
  mu_g_hat0 <- c(mu_g1_hat0, mu_g2_hat0, mu_g3_hat0)
  
  # 2. Get MLE under Hg1
  mu_g1_hat1 <- y_g1/C_1
  mu_g2_hat1 <- y_g2/C_2
  mu_g3_hat1 <- y_g3/C_3
  mu_g_hat1 <- c(mu_g1_hat1, mu_g2_hat1, mu_g3_hat1)
  
  # 3. Get test statistic
  logL <- function(y_set, C_set, mu_set){
    # if y_gi != 0, mu_gi != 0, for i=1,2,3
    if(mu_set[1] == 0){first <- 0}
    else{first <- y_set[1]*log(mu_set[1]) - C_set[1]*mu_set[1]}
    if(mu_set[2] == 0){second <- 0 }
    else{second <- y_set[2]*log(mu_set[2]) - C_set[2]*mu_set[2]}
    if(mu_set[3] == 0){third <- 0}
    else{third <- y_set[3]*log(mu_set[3]) - C_set[3]*mu_set[3]}
    return(first + second + third)
    #- log(factorial(y_set[1])) - log(factorial(y_set[2])) - log(factorial(y_set[3]))
    #+ y_set[1]log(C_set[1]) + y_set[2]log(C_set[2]) + y_set[3]log(C_set[3]) 
    # Since the factorial and ylogC terms will cancel even when y=0, but it can be infinity to make computation fail, so we omit them.
  }
  y_set = c(y_g1, y_g2, y_g3)
  C_set = c(C_1, C_2, C_3)
  LR <- -2 * (logL(y_set = y_set, C_set = C_set, mu_set = mu_g_hat0) - logL(y_set = y_set, C_set = C_set, mu_set = mu_g_hat1))
  pval <- pchisq(LR, df = 1, lower.tail = FALSE)
  # return(list(pval=pval, mu_g_hat0 = mu_g1_hat0, mu_g_hat1 = mu_g_hat1))
  return(pval)
}

Hetero_NB_LRT <- function(y_g1, y_g2, y_g3, C_1 = 1, C_2 = 2, C_3 = 1, dispersion = 1){
  # This function conducts NB LRT test for heterosis in a single family and for a single gene with single replication per genotype given dispersion parameter
  # Hg0: mu_g1 + mu_g2 = 2 * mu_g3
  # @param y_g1: raw count for female
  # @param y_g2: raw count for male
  # @param y_g3: raw count for hybrid
  # @param C_1: normalization factor for female sample
  # @param C_2: normalization factor for male sample
  # @param C_3: normalization factor for hybrid sample
  # @param dispersion: estimation of dispersion parameter
  # @output pval: pvalue for NB LRT test with known dispersion
  
  # rate = 1/dispersion
  if(dispersion == 0){return(NA)}
  r <- 1/dispersion
  
  # 1. Get MLE under global
  mu_g_hat1 <- c(y_g1/C_1, y_g2/C_2, y_g3/C_3)
  
  # 2. Get MLE under null
  y_set <- c(y_g1, y_g2, y_g3)
  C_set <- c(C_1, C_2, C_3)
  if(sum(y_set) == 0){return(1)}else{
    equations <- function(mu_set) {
      eq1 <- y_set[1]/mu_set[1] + y_set[3]/(mu_set[1]+mu_set[2]) - C_1*(r + y_set[1])/(C_1*mu_set[1] + r) - C_3*(r + y_set[3])/(C_3*(mu_set[1] + mu_set[2]) + 2*r)
      eq2 <- y_set[2]/mu_set[2] + y_set[3]/(mu_set[1]+mu_set[2]) - C_2*(r + y_set[2])/(C_2*mu_set[2] + r) - C_3*(r + y_set[3])/(C_3*(mu_set[1] + mu_set[2]) + 2*r)
      return(c(eq1, eq2))
    }
    
    init_list <- list(c(ifelse(y_g1 == 0, 10, y_g1), ifelse(y_g2 == 0, 10, y_g2)),
                      c(y_g1/C_1, y_g2/C_2), c(y_g1/2, y_g2/2),
                      c(10, 10), c(100, 100), c(2, 2), c(100, 10), c(10, 100),
                      c(100, 2), c(2, 100))
    
    # function to check the validity of the solution
    check_solution <- function(sol) {
      all(sol$x >= 0) && all(sol$x <= max(c(y_g1, y_g2, y_g3)) * 100) && sol$termcd < 4
    }
    
    # loop through the starting points to find a valid solution
    sol_found <- FALSE
    for (inits in init_list) {
      solution <- nleqslv(x = inits, fn = equations, control = list(maxit = 500))
      if (check_solution(solution)) {
        mu_g_hat0 <- c(solution$x, mean(solution$x))
        sol_found <- TRUE
        break
      }
    }
    
    # if no solution is found with nleqslv, use multiroot
    if (!sol_found) {
      root <- multiroot(f = equations, start = c(ifelse(y_g1 == 0, 10, y_g1/C_1), ifelse(y_g2 == 0, 10, y_g2/C_2)))
      if (all(root$root >= 0) && all(root$root <= max(c(y_g1, y_g2, y_g3)) * 100)) {
        mu_g_hat0 <- c(root$root, mean(root$root))
      } else {
        return(NA) # Or handle the failure to find a solution as appropriate
      }
    }
    
    logL <- function(y_set, C_set, mu_set, r){
      # under H1, if y_gi!=0, mu_gi!=0, under H0, its very rare to get 0 as MLE
      if(mu_set[1] == 0){first <- -(r + y_set[1])*log(r + C_set[1]*mu_set[1])}
      else{first <- y_set[1]*log(mu_set[1]) - (r + y_set[1])*log(r + C_set[1]*mu_set[1])}
      if(mu_set[2] == 0){second <- -(r + y_set[2])*log(r + C_set[2]*mu_set[2])}
      else{second <- y_set[2]*log(mu_set[2]) - (r + y_set[2])*log(r + C_set[2]*mu_set[2])}
      if(mu_set[3] == 0){third <- -(r + y_set[3])*log(r + C_set[3]*mu_set[3])}
      else{third <- y_set[3]*log(mu_set[3]) - (r + y_set[3])*log(r + C_set[3]*mu_set[3])}
      return(first + second + third)
      # Since the choose terms and rlogr and ylogC terms will cancel even with y=0, but it can be greater than the computation can afford  to make computation fail, so we omit them right now
    }
    LR <- -2 * (logL(y_set=y_set, C_set = C_set, mu_set = mu_g_hat0, r=r) - logL(y_set=y_set, C_set = C_set, mu_set = mu_g_hat1, r=r))
    pval <- pchisq(LR, df = 1, lower.tail = FALSE)
    # return(list(pval=pval, mu_g_hat0 = mu_g1_hat0, mu_g_hat1 = mu_g_hat1))
    return(pval)
  }
}

Pois_3mean_LRT <- function(y_g1, y_g2, y_g3, C_1 = 1, C_2 = 1, C_3 = 1){
  # This function conducts Poisson LRT test for non-family in a single family and for a single gene with single replication per genotype
  # Hg0: mu_g1 = mu_g2 = mu_g3
  # @param y_g1: raw count for female
  # @param y_g2: raw count for male
  # @param y_g3: raw count for hybrid
  # @param C_1: normalization factor for female sample
  # @param C_2: normalization factor for male sample
  # @param C_3: normalization factor for hybrid sample
  
  y_set <- c(y_g1, y_g2, y_g3)
  C_set <- c(C_1, C_2, C_3)
  mu_g_hat0 <- c(sum(y_set)/sum(C_set), sum(y_set)/sum(C_set), sum(y_set)/sum(C_set))
  mu_g_hat1 <- c(y_g1/C_1, y_g2/C_2, y_g3/C_3)
  logL <- function(y_set, C_set, mu_set){
    # if y_gi != 0, mu_gi != 0, for i=1,2,3
    if(mu_set[1] == 0){first <- 0}
    else{first <- y_set[1]*log(mu_set[1]) - C_set[1]*mu_set[1]}
    if(mu_set[2] == 0){second <- 0 }
    else{second <- y_set[2]*log(mu_set[2]) - C_set[2]*mu_set[2]}
    if(mu_set[3] == 0){third <- 0}
    else{third <- y_set[3]*log(mu_set[3]) - C_set[3]*mu_set[3]}
    return(first + second + third)
    #- log(factorial(y_set[1])) - log(factorial(y_set[2])) - log(factorial(y_set[3]))
    #+ y_set[1]log(C_set[1]) + y_set[2]log(C_set[2]) + y_set[3]log(C_set[3]) 
    # Since the factorial and ylogC terms will cancel even when y=0, but it can be infinity to make computation fail, so we omit them.
  }
  LR <- -2 * (logL(y_set = y_set, C_set = C_set, mu_set = mu_g_hat0) - logL(y_set = y_set, C_set = C_set, mu_set = mu_g_hat1))
  pval <- pchisq(LR, df = 2, lower.tail = FALSE)
  # return(list(pval=pval, mu_g_hat0 = mu_g1_hat0, mu_g_hat1 = mu_g_hat1))
  return(pval)
}

Disp_est_null_family <- function(counts, group = NULL){
  # This function estimate the tagwise dispersion parameter for multiple genes with multiple family
  # @param counts: with genes in rows, and three*n_families columns, female, male and hybrid which passed non-family LRT
  # @param group: a vector indicating the group information of three*n_families columns
  if(is.null(group)){
    y <- DGEList(counts = counts)
    y <- calcNormFactors(y, method = "RLE")
    y <- estimateCommonDisp(y)
    y <- estimateTrendedDisp(y)
    y <- estimateTagwiseDisp(y)
  }else{
    y <- DGEList(counts = counts, group = group) 
    y <- calcNormFactors(y, method = "RLE")
    design <- model.matrix(~ group)
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTrendedDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
  }
  return(y$tagwise.dispersion)
}

twoStageLRT <- function(hybrid_mat, female_mat, male_mat, disp_mat = NULL,
                        normalization = FALSE, norm_factors = NULL,
                        n_fam_thres = 20, gene_groups_no = 3){
  # This function apply two stage LRT to data with formats: three count matrices with genes in rows, samples in columns
  # @param hybrid_mat: hybrids expression count with genes in rows, samples in columns.
  # @param female_mat: female expression count with genes in rows, samples in columns.
  # @param male_mat: male expression count with genes in rows, samples in columns.
  # @param disp_mat: If "Pois", use Pois LRT and only output Pois results
  #                  if NULL, use dispersion estimation, else, use provided disp_mat for NB test, but output will have both NBLRT and PoisLRT results
  # @param normalization: whether or not need normalization for the raw data matrices, 
  #                       if FALSE, treat the normalization factors all 1.
  # @param norm_factors: if normalization = TRUE and norm_factors = NULL, use DESeq normalization factors,
  #                      if normalization = TRUE and norm_factors is a list contains three vectors 
  #                      female, male and hybrid normalization factors, use the provided norm_factors in tests
  # @param n_fam_thres: genes with at least n_fam_thres null families will be used for dispersion estimation, 
  #                     otherwise, will use loess prediction to get dispersion, default is 20
  # @param gene_groups_no: the number of clusters for genes, default is 3.
  # @output: if disp_mat = NULL, return Poisson LRT for heterosis pvalue and qvalue from BH adjustment
  #          if disp_mat != NULL, return both NB and Poisson LRT for heterosis pvalue and qvalue and genewise dispersion.
  #          and the estimation of mu_gi, which is the female/male/hybrid_mat_normed
  
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
    
  }else if(normalization & (!is.null(norm_factors))){
    cat("Use provided normalization factors, no estimation. \n")
    female_mat_normed <- sweep(female_mat, 2, norm_factors$female_norm_factor, "/")
    male_mat_normed <- sweep(male_mat, 2, norm_factors$male_norm_factor, "/")
    hybrid_mat_normed <- sweep(hybrid_mat, 2, norm_factors$hybrid_norm_factor, "/")
  }else{ # when normalization is FALSE, all norm_factors will be set to be 1
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
    print(dim(null_disp_ind_sub))
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
  }else if(! "character" %in% class(disp_mat)){ # use provided dispersion parameters
    dispersion_mat <- disp_mat
  }else{ # use Poisson LRT for MPH
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

Heter_pip_2sLRT <- function(n_genes, n_families,
                            mu_female, mu_male, effect_factor,
                            female_norm_factor, male_norm_factor, hybrid_norm_factor, true_norm = c(TRUE, FALSE),
                            MP_DISP_func, true_disp = c("None", TRUE, FALSE), 
                            n_fam_thres = 20, gene_groups_no = 3, repetition = 150, 
                            plot = TRUE, filename = NULL){
  # This function simulates n_families with each family includes three counts (i=1-female, 2-male, 3-hybrid) from Negative Binomial distribution for n_genes. Then two stage LRT is applied to check heterosis gene-family combination 
  # G genes (specify) with normalization factors C_i, mean mu_{gi} and dispersion parameter phi_g for each family.
  # Hg0: mu_{g1} + mu_{g2} = 2mu_{g3} gene g is non heterosis in this family.
  # Hg1: otherwise. The effect size under alternative is specified by effect factor - how many times more mu_g3 is of (mu_g1 + mu_g2) / 2, 0 = Hg0
  # @param n_genes: number of genes.
  # @param n_families: number of families.
  # @param female_norm_factor: a vector with normalization factors for n_families females.
  # @param hybrid_norm_factor: a vector with normalization factors for n_families males.
  # @param male_norm_factor: a vector with normalization factors for n_families hybrids.
  #         if norm_factor == NULL, norm_factor will be all 1.
  #         if provided and are not all 1, DESeq normalization factors will be used in 2sLRT test.
  #         if provided with length = 2, norm_factor will be simulated from Unif(norm_factor[1], norm_factor[2]), DESeq will also be used.
  #         if provided with length = 1, all samples share the same normalization factor, DESeq will not be used.
  # @param MP_DISP_func: character name of function describing the relationship between mid-parent and dispersion parameters.
  # @param mu_female: a matrix containing mean expression count of female samples for all genes and all families with nrow = n_genes, ncol = n_families.
  # @param mu_male: a matrix containing mean count of male samples for n_genes and n_families.
  #         if a vector with length = 2, mean mu (not Cmu) for all genes and all families in male sample will be simulated from Gamma(mu_male[1], mu_male[2]) for shape and rate, thus the mean of gamma distribution is shape/rate.
  #         if a vector with length = 1, all male samples share the same mean count.
  # @param effect_factor: a matrix containing how many times more mu_{g3} is of 
  #         (mu_{g1} + mu_{g2})/2, for null cases, set it to be 0. 
  #         if a vector with length = 1, all families share the same effect size.
  # @param true_disp: "None" for Pois LRT, TRUE for NB LRT with true disp, FALSE for NB LRT with estimated disp.
  # @param true_norm: TRUE for using true C_i in tests, FALSE for using DESeq normalization factors in tests
  # @param n_fam_thres: genes with at least n_fam_thres all-null families will be used for dispersion estimation, otherwise, will use loess prediction to get dispersion, default is 20
  # @param gene_groups_no: the number of clusters for genes, default is 3.
  # @param repetition: the number of replications for simulation
  # @param filename: if NULL, no save of any plots, if specified, means on slurm running, will save plots.
  # @output: a list including pvalue matrix and qvalue matrix (computed using BH method).
  
  
  
  # mu_female: n_genes * n_families
  if(length(mu_female) == 1){
    cat("All genes and families have the same mean female count.\n")
    mu_female <- matrix(mu_female, nrow = n_genes, ncol = n_families)
    }else if(length(mu_female) == 2){
    cat("Use simulated mean count from gamma with shape and rate for female. \n")
    mu_female <- matrix(rgamma(n_genes * n_families, shape = mu_female[1], rate = mu_female[2]), nrow = n_genes, byrow = T)
    }else if(nrow(mu_female) == n_genes & ncol(mu_female) == n_families){
    cat("Use provided mean count matrix for female. \n")
        }else{cat("Error: mu_female must be a matrix or a vector of length 2. \n")}
  
  # mu_male: n_genes * n_families
  if(length(mu_male) == 1){
    cat("All genes and families have the same mean male count. \n")
    mu_male <- matrix(mu_male, nrow = n_genes, ncol = n_families)
    }else if(length(mu_male) == 2){
        cat("Use simulated mean count from gamma with shape and rate for male. \n")
        mu_male <- matrix(rgamma(n_genes * n_families, shape = mu_male[1], rate = mu_male[2]), nrow = n_genes, byrow = T)
        }else if(nrow(mu_male) == n_genes & ncol(mu_male) == n_families){
            cat("Use provided mean count matrix for male. \n")
        }else{
                cat("Error: mu_male must be a matrix or a vector of length 2. \n")
            }
  
  # effect_factor: n_genes * n_families
  if(length(effect_factor) == 1){
    effect_factor <- matrix(effect_factor, nrow = n_genes, ncol = n_families)
    }else if(nrow(effect_factor) == n_genes & ncol(effect_factor) == n_families){
    cat("Use provided effect factor matrix for hybrid. \n")
    }else{
            cat("Error: effect factor must have length either 1 or n_genes. \n")
        }
  
  # mu_hybrid: n_genes * n_families
  mu_hybrid <- (effect_factor + 1)*((mu_female + mu_male)/2)
  
  # female_norm_factors:  vector with size = n_families
  # male_norm_factors:  vector with size = n_families
  # hybrid_norm_factors:  vector with size = n_families
  norm_generator <- function(norm_factor){
    if(is.null(norm_factor)){
      cat("Use all 1 normalization factors for all samples. \n")
      norm_factor <- rep(1, n_families)
      }else if(length(norm_factor) == 2){
      cat("Use simulated norm factors from uniform. \n")
      norm_factor <- c(runif(n_families, min = norm_factor[1], max = norm_factor[2]))
      }else if(length(norm_factor) == 1){
      norm_factor <- rep(norm_factor, n_families)
      }else{
          cat("Error: norm_factors must have either 2 or 1 elements. \n")
          }
    return(norm_factor)
  }
  female_norm_factor <- norm_generator(female_norm_factor)
  male_norm_factor <- norm_generator(male_norm_factor)
  hybrid_norm_factor <- norm_generator(hybrid_norm_factor)
  if(mean(female_norm_factor==1) == 1 & mean(male_norm_factor==1) == 1 & 
     mean(hybrid_norm_factor == 1) == 1){
      normalization <- FALSE
      }else{normalization <- TRUE} 
  # normalization indicates whether we will have estimated norm factors, but whether use estimation is decided by true_norm
  
  # lambda_gi = C_imu_gi
  female_lambda <- sweep(mu_female, 2, female_norm_factor, `*`)
  male_lambda <- sweep(mu_male, 2, male_norm_factor, `*`)
  hybrid_lambda <- sweep(mu_hybrid, 2, hybrid_norm_factor, `*`)
  
  # dispersion: n_genes * n_families
  mp <- (mu_female + mu_male)/2
  mp_seeds <- apply(mp, 1, median)
  if(MP_DISP_func == "low"){ # use the function to simulate
      mp_seeding_disps <- sapply(mp_seeds, inverse_func)
  }else if(MP_DISP_func == "high"){
      mp_seeding_disps <- get_sim_dispersion(mp_seeds, sim_params_med, scale_noise=1)
  }else if(MP_DISP_func == "extreme"){
      mp_seeding_disps <- get_sim_dispersion(mp_seeds, sim_params_med, scale_noise=1.5)
  }else{
      stop("Invalid MP_DISP_func!")
  }
  
  # gene wise dispersion
  dispersion <- matrix(mp_seeding_disps, nrow = n_genes, ncol = n_families)
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
    if(true_disp == "None"){
        disp_mat <- "Pois" # only use Poisson LRT for heterosis
    } else if(true_disp){
        disp_mat <- dispersion # use NB LRT with true disp and Pois LRT for heterosis
        } else{
            disp_mat <- NULL
            } # use NB LRT with estimated disp and Pois LRT for heterosis
    
    # determine use what type of normalization factors in test
    if(true_norm){norm_factors <- list(female_norm_factor=female_norm_factor, # use true normalization factors in test if normalization = TRUE
                                       male_norm_factor=male_norm_factor,
                                       hybrid_norm_factor=hybrid_norm_factor)
    } else{norm_factors <- NULL} # use DESeq estimated norm factors if normalization = TRUE
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

      # dispersion bias and RMSE generation
      diff_matrix <- sweep(DISP_est, 2, mp_seeding_disps, FUN = "-")
      # BIAS: Average error across all genes for each run
      # (We use rowMeans because each row is a simulation run)
      bias_per_run <- rowMeans(diff_matrix, na.rm = TRUE)
      # RMSE: Square root of the mean squared error across all genes for each run
      rmse_per_run <- sqrt(rowMeans(diff_matrix^2, na.rm = TRUE))
      n_runs <- length(bias_per_run)
      # Average Bias and Standard Error (across the repetitions)
      final_bias_mean <- mean(bias_per_run)
      final_bias_se   <- sd(bias_per_run) / sqrt(n_runs)
      # Average RMSE and Standard Error
      final_rmse_mean <- mean(rmse_per_run)
      final_rmse_se   <- sd(rmse_per_run) / sqrt(n_runs)
      cat(sprintf("Bias: %.5f (SE: %.5f)\n", final_bias_mean, final_bias_se))
      cat(sprintf("RMSE: %.5f (SE: %.5f)\n", final_rmse_mean, final_rmse_se))
      
      # --- ADDITIONAL median-based summaries (no changes above) ---
      
      # 1) Median bias across genes within each run (robust per-run bias)
      bias_per_run_med_genes <- apply(diff_matrix, 1, median, na.rm = TRUE)
      
      # 3) Median-of-medians: median across runs of median-over-genes
      final_bias_med_genes_med_runs <- median(bias_per_run_med_genes, na.rm = TRUE)
      final_bias_med_genes_mean_runs <- mean(bias_per_run_med_genes, na.rm = TRUE)
      
      
      cat(sprintf("Median-of-medians bias (median over genes, then median over runs): %.5f\n",
                  final_bias_med_genes_med_runs))
      cat(sprintf("Mean-of-medians bias (median over genes, then mean over runs): %.5f\n",
                  final_bias_med_genes_mean_runs))
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
    }else{
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
      
      #cat(paste("\n Average TPR_max with empirical FDR <= 0.05 is, NB=",
      #          round(at_actualFDR(qval_NB),4),
      #          ", Pois=",  round(at_actualFDR(qval_Pois),4), "\n", sep =""))
      
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
      
    }else{
      p_roc_avg <- NULL
      p_FDR_avg <- NULL
      p_TPR_avg <- NULL
      AUC_avg <- NA
    }
    if(is.null(qq0)){
        finalp <- plot_grid(p_disp, p_disp_used, p_disp_scatter, p_single_disp_vis, 
                                         p_roc_avg, p_FDR_avg, p_TPR_avg, nrow = 2, ncol = 4)
    }else if(is.null(p_roc_avg)){
        finalp <- plot_grid(p_disp, p_disp_used, p_disp_scatter, p_single_disp_vis, qq0, nrow = 2)
        }else{
      finalp <- plot_grid(p_disp, p_disp_used, p_disp_scatter, p_single_disp_vis, 
                          qq0, p_roc_avg, p_FDR_avg, p_TPR_avg, nrow = 2, ncol = 4)
      paperp <- plot_grid(qq0, p_roc_avg, p_FDR_avg, p_TPR_avg, nrow = 1)
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

inverse_func <- function(mp){
  disp <- (mp)^(-1/6) - 0.3 + rnorm(1, 0, sd = 0.002)
  return(disp)
}

# -------------- Running  for equal parental mean -------------
# load("/rafalab/yunhui/Others/MoreDispersion.RData") # for more dispersion patterns like high and medium, borrowed from edgeR manual
# if(pattern_arg == "dense"){
#   # dense heterosis pattern
#   mu_female <- "matrix(rgamma(case$n_genes*case$n_families, shape = 30, rate=runif(case$n_genes, 0.05, 0.2)),nrow=case$n_genes, ncol=case$n_families)"
#   effect_factor <- "rbind(matrix(ifelse(rbinom(case$n_genes/4*case$n_families, 1, 0.05)==1, case$sig_pos, 0), nrow = case$n_genes/4, ncol = case$n_families), matrix(ifelse(rbinom(case$n_genes/4*case$n_families, 1, 0.05)==1, case$sig_neg, 0),nrow = case$n_genes/4, ncol = case$n_families), matrix(ifelse(rbinom(case$n_genes/4*case$n_families, 1, 0.9)==1, case$sig_pos, 0),nrow = case$n_genes/4, ncol = case$n_families), matrix(ifelse(rbinom(case$n_genes/4*case$n_families, 1, 0.9)==1, case$sig_neg, 0),nrow = case$n_genes/4, ncol = case$n_families))"
# }else{
#   # sparse heterosis pattern
#   load("/rafalab/yunhui/Others/Normfac_Prop_sig_hyb_genewise_20284.RData")
#   mu_female <- "matrix(rgamma(case$n_genes*case$n_families, shape = 30, rate=runif(case$n_genes, 0.05, 0.2)), nrow=case$n_genes, ncol=case$n_families)"
#   effect_factor <- "rbind(matrix(ifelse(rbinom(case$n_genes/2*case$n_families, 1, sample(prop_sig_hybrids_per_gene, size = case$n_genes/2, replace=TRUE)) == 1, case$sig_pos, 0), nrow = case$n_genes/2, ncol = case$n_families), 
#   matrix(ifelse(rbinom(case$n_genes/2*case$n_families, 1, sample(prop_sig_hybrids_per_gene, size = case$n_genes/2, replace=TRUE)) == 1, case$sig_neg, 0), nrow = case$n_genes/2, ncol = case$n_families))"
# }
# 
# 
# 
# case <- list(n_genes = ngenes_arg, n_families = nfamilies_arg,
#              sig_pos = 3, sig_neg = -3/4,
#              mu_female = mu_female,
#              mu_male = "case$mu_female",
#              effect_factor = effect_factor,
#              female_norm_factor = c(0.8, 1.3), male_norm_factor = c(0.8, 1.3), hybrid_norm_factor = c(0.8, 1.3),
#              MP_DISP_func = "inverse_func", true_disp = FALSE, true_norm = FALSE,
#              n_fam_thres = N0_arg, gene_groups_no = K_arg, repetition = 10, plot = TRUE, filename = output_filename_arg)
# print(case)
# if(is.character(case$mu_female)){case$mu_female <- eval(parse(text = case$mu_female))}
# if(is.character(case$mu_male)){case$mu_male <- eval(parse(text = case$mu_male))}
# if(is.character(case$effect_factor)){case$effect_factor <- eval(parse(text = case$effect_factor))}
# res <- Heter_pip_2sLRT(n_genes = case$n_genes, n_families = case$n_families,
#                        mu_female = case$mu_female, mu_male = case$mu_male,
#                        effect_factor = case$effect_factor,
#                        female_norm_factor = case$female_norm_factor,
#                        male_norm_factor = case$male_norm_factor,
#                        hybrid_norm_factor = case$hybrid_norm_factor,
#                        true_norm = case$true_norm,
#                        MP_DISP_func = case$MP_DISP_func, true_disp = case$true_disp,
#                        n_fam_thres = case$n_fam_thres,
#                        gene_groups_no = case$gene_groups_no,
#                        repetition = case$repetition, plot = case$plot,
#                        filename = case$filename)
# 


################## For unequal parental mu for the second half of the genes ##################
load("/rafalab/yunhui/Others/MoreDispersion.RData") # for more dispersion patterns like high and medium, borrowed from edgeR manual
if(pattern_arg == "dense"){
    # dense heterosis pattern
    effect_factor <- "rbind(matrix(ifelse(rbinom(case$n_genes/4*case$n_families, 1, 0.05)==1, case$sig_pos, 0),
                                   nrow = case$n_genes/4, ncol = case$n_families), 
                            matrix(ifelse(rbinom(case$n_genes/4*case$n_families, 1, 0.05)==1, case$sig_neg, 0),
                                   nrow = case$n_genes/4, ncol = case$n_families), 
                            matrix(ifelse(rbinom(case$n_genes/4*case$n_families, 1, 0.9)==1, case$sig_pos, 0),
                                   nrow = case$n_genes/4, ncol = case$n_families), 
                            matrix(ifelse(rbinom(case$n_genes/4*case$n_families, 1, 0.9)==1, case$sig_neg, 0),
                                   nrow = case$n_genes/4, ncol = case$n_families))"
}else{
    # sparse heterosis pattern
    load("/rafalab/yunhui/Others/Normfac_Prop_sig_hyb_genewise_20284.RData")
    effect_factor <- "rbind(matrix(ifelse(rbinom(case$n_genes/2*case$n_families, 1, sample(prop_sig_hybrids_per_gene, size = case$n_genes/2, replace=TRUE)) == 1, case$sig_pos, 0), 
                                   nrow = case$n_genes/2, ncol = case$n_families), 
                                   matrix(ifelse(rbinom(case$n_genes/2*case$n_families, 1, sample(prop_sig_hybrids_per_gene, size = case$n_genes/2, replace=TRUE)) == 1, case$sig_neg, 0),
                                   nrow = case$n_genes/2, ncol = case$n_families))"
    

}
rate_genes <- "runif(case$n_genes, 0.05, 0.2)"
mu_female <- "matrix(rgamma(case$n_genes*case$n_families, 
                                shape = 30, 
                                rate=case$rate_genes),
                         nrow=case$n_genes, ncol=case$n_families)"
mu_male <- "rbind(case$mu_female[1:(case$n_genes/2),],
                      matrix(rgamma(case$n_genes*case$n_families/2, 
                                    shape = 30, 
                                    rate=case$rate_genes[(case$n_genes/2+1):case$n_genes]), nrow=case$n_genes/2, ncol=case$n_families))"

case <- list(n_genes = ngenes_arg, n_families = nfamilies_arg,
             sig_pos = 3, sig_neg = -3/4, # 2, -2/3 for weak, 3, -3/4 for strong
             rate_genes = rate_genes,
             mu_female = mu_female,
             mu_male = mu_male, 
             effect_factor = effect_factor,
             female_norm_factor = c(0.8, 1.3), male_norm_factor = c(0.8, 1.3), hybrid_norm_factor = c(0.8, 1.3),
             MP_DISP_func = disp_arg, true_disp = FALSE, true_norm = FALSE,
             n_fam_thres = N0_arg, gene_groups_no = K_arg, repetition = 10, plot = TRUE, filename = output_filename_arg)

print(case)
if(is.character(case$rate_genes)){case$rate_genes <- eval(parse(text = case$rate_genes))}
if(is.character(case$mu_female)){case$mu_female <- eval(parse(text = case$mu_female))}
if(is.character(case$mu_male)){case$mu_male <- eval(parse(text = case$mu_male))}
if(is.character(case$effect_factor)){case$effect_factor <- eval(parse(text = case$effect_factor))}
res <- Heter_pip_2sLRT(n_genes = case$n_genes, n_families = case$n_families,
                       mu_female = case$mu_female, mu_male = case$mu_male,
                       effect_factor = case$effect_factor,
                       female_norm_factor = case$female_norm_factor,
                       male_norm_factor = case$male_norm_factor,
                       hybrid_norm_factor = case$hybrid_norm_factor,
                       true_norm = case$true_norm,
                       MP_DISP_func = case$MP_DISP_func, true_disp = case$true_disp,
                       n_fam_thres = case$n_fam_thres,
                       gene_groups_no = case$gene_groups_no,
                       repetition = case$repetition, plot = case$plot,
                       filename = case$filename)


# ------  Code for generating the dispersions from edgeR manual

#' library(edgeR)
#' library(tweeDEseqCountData)
#' 
#' # =========================================================
#' # 1. PREPARE DATASETS (Arabidopsis)
#' # =========================================================
#' 
#' # --- ARABIDOPSIS (Medium Dispersion) ---
#' url <- "https://bioinf.wehi.edu.au/edgeR/UserGuideData/arab.rds"
#' if(!file.exists("arab.rds")) download.file(url, "arab.rds", mode="wb")
#' arab <- readRDS("arab.rds")
#' Treat <- factor(substring(colnames(arab),1,4))
#' Treat <- relevel(Treat, ref="mock")
#' Time <- factor(substring(colnames(arab),5,5))
#' y_med <- DGEList(counts=arab, group=Treat)
#' keep <- filterByExpr(y_med)
#' y_med <- y_med[keep, , keep.lib.sizes=FALSE]
#' y_med <- normLibSizes(y_med)
#' design_med <- model.matrix(~Time+Treat)
#' y_med <- estimateDisp(y_med, design_med, robust=TRUE)
#' 
#' 
#' # =========================================================
#' # 2. TRAINING FUNCTION (Internally uses Log Scale for Stability)
#' # =========================================================
#' get_simulation_parameters <- function(dge_object) {
#'     # We use LogCPM because it creates a stable, realistic curve
#'     means_log <- dge_object$AveLogCPM
#'     disps_log <- log(dge_object$tagwise.dispersion)
#'     
#'     # Calculate average library size (needed to convert your Raw MP later)
#'     avg_lib_size <- mean(dge_object$samples$lib.size)
#'     
#'     # Fit the trend
#'     df <- data.frame(mean = means_log, log_disp = disps_log)
#'     fit <- loess(log_disp ~ mean, data=df, span=0.3)
#'     
#'     # Calculate residuals (the noise)
#'     resids <- residuals(fit)
#'     
#'     return(list(fit=fit, residuals=resids, avg_lib_size=avg_lib_size))
#' }
#' 
#' sim_params_med  <- get_simulation_parameters(y_med)
#' 
#' 
#' # =========================================================
#' # 3. GENERATOR FUNCTION (Accepts YOUR Raw Means)
#' # =========================================================
#' 
#' #' Generate dispersions for simulated RAW Mid-Parent Means
#' #' 
#' #' @param raw_mp_means Vector of simulated RAW means (e.g., 10, 500, 2000...)
#' #' @param sim_params The object created above (high or med)
#' #' @return Vector of dispersion parameters
# get_sim_dispersion <- function(raw_mp_means, sim_params, scale_noise=1) {
# 
#     # 1. Convert Raw Means to LogCPM
#     log_cpm_means <- log2((raw_mp_means + 0.25) / sim_params$avg_lib_size * 1e6)
# 
#     # 2. Predict Trend
#     pred_log_disp <- predict(sim_params$fit, newdata=data.frame(mean=log_cpm_means))
# 
#     # 3. Handle Extrapolation
#     min_real <- min(sim_params$fit$x); max_real <- max(sim_params$fit$x)
#     if(any(is.na(pred_log_disp))) {
#         low_val <- predict(sim_params$fit, data.frame(mean=min_real))
#         high_val <- predict(sim_params$fit, data.frame(mean=max_real))
#         pred_log_disp[is.na(pred_log_disp) & log_cpm_means < min_real] <- low_val
#         pred_log_disp[is.na(pred_log_disp) & log_cpm_means > max_real] <- high_val
#     }
# 
#     # 4. Generate Initial Dispersions
#     n <- length(raw_mp_means)
#     random_resid <- sample(sim_params$residuals, size=n, replace=TRUE)
# 
#     # Calculate in Log scale first
#     current_log_val <- pred_log_disp + (random_resid * scale_noise)
# 
#     # 5. Resampling Loop (Ensure Disp <= 1)
#     bad_idx <- which(current_log_val > log(0.75))
#     iter <- 0
#     max_iter <- 50 # Prevent infinite loops if trend is extremely high
# 
#     while(length(bad_idx) > 0 && iter < max_iter) {
#         # Resample residuals ONLY for the problematic genes
#         new_resids <- sample(sim_params$residuals, size=length(bad_idx), replace=TRUE)
# 
#         # Recalculate
#         current_log_val[bad_idx] <- pred_log_disp[bad_idx] + (new_resids * scale_noise)
# 
#         # Re-check
#         bad_idx <- which(current_log_val > log(0.75))
#         iter <- iter + 1
#     }
# 
#     # Safety valve: If we hit max_iter and some are still > 1, clamp them to 1.0
#     if(length(bad_idx) > 0) {
#         current_log_val[bad_idx] <- 0 # exp(0) == 1
#     }
# 
#     return(exp(current_log_val))
# }
#' # =========================================================
#' # 4. EXAMPLE USAGE IN YOUR SIMULATION
#' # =========================================================
# case$n_genes <- 10000
# case$n_fam_thres <- 150
# rate_genes <- runif(case$n_genes, 0.05, 0.2)
# mu_female <- matrix(rgamma(case$n_genes*case$n_families,
#                                 shape = 30,
#                                 rate=case$rate_genes),
#                          nrow=case$n_genes, ncol=case$n_families)
# mu_male <- rbind(case$mu_female[1:(case$n_genes/2),],
#                       matrix(rgamma(case$n_genes*case$n_families/2,
#                                     shape = 30,
#                                     rate=case$rate_genes[(case$n_genes/2+1):case$n_genes]), nrow=case$n_genes/2, ncol=case$n_families))
# 
# mp <- (mu_female + mu_male)/2
# mp_seeds <- apply(mp, 1, median)
# df <- data.frame(MedianMP=mp_seeds,
#                  Disp_low=sapply(mp_seeds, inverse_func),
#                  Disp_high=get_sim_dispersion(mp_seeds, sim_params_med, scale_noise=1),
#                  Disp_extreme = get_sim_dispersion(mp_seeds, sim_params_med, scale_noise=1.5))
# df <- tidyr::pivot_longer(df, cols = c("Disp_low","Disp_high","Disp_extreme"), names_to="Disp_level", values_to = "Disp")
# df$Disp_level <- factor(df$Disp_level, levels = c("Disp_low","Disp_high","Disp_extreme"), labels = c("Low","High","Extreme"))
# p <- ggplot(df)+geom_point(aes(x=MedianMP, y= Disp, color = Disp_level))+facet_wrap(vars(Disp_level))+theme_bw()+
#     theme(legend.position = "bottom")+labs(x="Median mid-parental mean", y = "Dispersion parameters", color = "Variability Level")
# ggsave(filename = "DispersionLevelExt.pdf",plot = p, width = 10, height = 4, units = "in")

