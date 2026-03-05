library(ggplot2)
library(pROC)
library(limma)
library(MASS)
library(ggpubr)
library(nleqslv)
library(rootSolve, lib.loc="/rafalab/yunhui/Rpackages/") # Ensure lib.loc is correct for your system if needed
library(pracma)
library(edgeR)
library(DESeq2)

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


inverse_func <- function(mp){
  disp <- (mp)^(-1/6) - 0.3 + rnorm(1, 0, sd = 0.002)
  return(disp)
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
  
  # # Poisson LRT for heterosis
  # pval_Pois <- matrix(mapply(Hetero_Pois_LRT, female_mat, male_mat, hybrid_mat,
  #                            female_norm_factor_mat, male_norm_factor_mat, hybrid_norm_factor_mat), 
  #                     nrow = n_gene, ncol = n_family)
  # rownames(pval_Pois) <- rownames(hybrid_mat)
  # colnames(pval_Pois) <- colnames(hybrid_mat)
  # qval_Pois <- apply(pval_Pois, 2, function(x){p.adjust(x, method = "fdr")})
  pval_Pois <- qval_Pois <- NULL
  
  
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

# average FDR and TPR curve of qvalue with nominal alpha in (0,0.2)
alpha_func <- function(qvalues, effect_factor){
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

at_actualFDR <- function(qvalues, effect_factor, actual_FDR_level = 0.05){
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

# average and SE AUC and partial AUC only for FPR in 0,0.2 and mean ROC curve
mean_ROC_func <- function(qvalues, effect_factor, direction = c(">", "<")){
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

# --- NEW: Replicated Newton-Raphson Logic ---

logLik_NB <- function(y, C, mu, r) { sum(y * log(C * mu) - (y + r) * log(C * mu + r)) }
score_NB  <- function(y, C, mu, r) { sum( (y / mu) - ((y + r) * C) / (C * mu + r) ) }
hessian_NB<- function(y, C, mu, r) { sum( -(y / mu^2) + ((y + r) * C^2) / (C * mu + r)^2 ) }

solve_H1_1D <- function(y, C, r) {
  mu <- sum(y) / sum(C)
  if(is.nan(mu) || mu == 0) return(1e-4)
  for(k in 1:25) {
    S <- score_NB(y, C, mu, r); H <- hessian_NB(y, C, mu, r)
    step <- S / H
    mu_new <- mu - step
    iter_half <- 0
    while(mu_new <= 0 && iter_half < 10) { step <- step * 0.5; mu_new <- mu - step; iter_half <- iter_half + 1 }
    if(mu_new <= 0) mu_new <- 1e-5
    if(abs(mu_new - mu) < 1e-6 * mu) return(mu_new)
    mu <- mu_new
  }
  return(mu)
}

test_Family_LRT <- function(y_list, C_list, phi) {
  r <- 1/phi
  mu_hat <- numeric(3); logL_H1 <- 0
  for(i in 1:3) {
    mu_hat[i] <- solve_H1_1D(y_list[[i]], C_list[[i]], r)
    logL_H1 <- logL_H1 + logLik_NB(y_list[[i]], C_list[[i]], mu_hat[i], r)
  }
  
  theta <- c(mu_hat[1], mu_hat[2])
  for(iter in 1:25) {
    m1 <- theta[1]; m2 <- theta[2]; m3 <- (m1 + m2)/2
    d1 <- c(score_NB(y_list[[1]],C_list[[1]],m1,r), hessian_NB(y_list[[1]],C_list[[1]],m1,r))
    d2 <- c(score_NB(y_list[[2]],C_list[[2]],m2,r), hessian_NB(y_list[[2]],C_list[[2]],m2,r))
    d3 <- c(score_NB(y_list[[3]],C_list[[3]],m3,r), hessian_NB(y_list[[3]],C_list[[3]],m3,r))
    
    grad <- c(d1[1] + 0.5*d3[1], d2[1] + 0.5*d3[1])
    H11 <- d1[2] + 0.25*d3[2]; H22 <- d2[2] + 0.25*d3[2]; H12 <- 0.25*d3[2]
    
    det <- H11*H22 - H12*H12
    if(abs(det) < 1e-10) break
    delta1 <- (H22*grad[1] - H12*grad[2])/det
    delta2 <- (-H12*grad[1] + H11*grad[2])/det
    
    theta_new <- theta - c(delta1, delta2)
    step_size <- 1
    while(any(theta_new <= 0) && step_size > 1e-4) {
      step_size <- step_size * 0.5
      theta_new <- theta - (c(delta1, delta2) * step_size)
    }
    if(max(abs(theta_new - theta)) < 1e-5) { theta <- theta_new; break }
    theta <- theta_new
  }
  logL_H0 <- logLik_NB(y_list[[1]],C_list[[1]],theta[1],r) + 
    logLik_NB(y_list[[2]],C_list[[2]],theta[2],r) + 
    logLik_NB(y_list[[3]],C_list[[3]],(theta[1]+theta[2])/2,r)
  
  LRT <- 2 * (logL_H1 - logL_H0)
  if(LRT < 0 || is.na(LRT)) LRT <- 0
  return(pchisq(LRT, df=1, lower.tail=FALSE))
}
# ==============================================================================
# 5. MAIN SIMULATION FUNCTION (UNREPLICATED vs REPLICATED)
# ==============================================================================

Heter_pip_Comparison <- function(n_genes, n_families, n_reps,
                                 mu_female, mu_male, effect_factor,
                                 female_norm_factor, male_norm_factor, hybrid_norm_factor,
                                 MP_DISP_func, n_fam_thres = 20, gene_groups_no = 3, 
                                 repetition = 50, output_file = NULL) {
  
  # --- 1. SETUP PARAMETERS (Your Original Code) ---
  if(length(mu_female) == 1){
    cat("All genes and families have the same mean female count.\n")
    mu_female <- matrix(mu_female, nrow = n_genes, ncol = n_families)
  }else if(length(mu_female) == 2){
    cat("Use simulated mean count from gamma with shape and rate for female. \n")
    mu_female <- matrix(rgamma(n_genes * n_families, shape = mu_female[1], rate = mu_female[2]), nrow = n_genes, byrow = T)
  }else if(nrow(mu_female) == n_genes & ncol(mu_female) == n_families){
    cat("Use provided mean count matrix for female. \n")
  }else{cat("Error: mu_female must be a matrix or a vector of length 2. \n")}
  
  if(length(mu_male) == 1){
    cat("All genes and families have the same mean male count. \n")
    mu_male <- matrix(mu_male, nrow = n_genes, ncol = n_families)
  }else if(length(mu_male) == 2){
    cat("Use simulated mean count from gamma with shape and rate for male. \n")
    mu_male <- matrix(rgamma(n_genes * n_families, shape = mu_male[1], rate = mu_male[2]), nrow = n_genes, byrow = T)
  }else if(nrow(mu_male) == n_genes & ncol(mu_male) == n_families){
    cat("Use provided mean count matrix for male. \n")
  }else{cat("Error: mu_male must be a matrix or a vector of length 2. \n")}
  
  if(length(effect_factor) == 1){
    effect_factor <- matrix(effect_factor, nrow = n_genes, ncol = n_families)
  }else if(nrow(effect_factor) == n_genes & ncol(effect_factor) == n_families){
    cat("Use provided effect factor matrix for hybrid. \n")
  }else{cat("Error: effect factor must have length either 1 or n_genes. \n")}
  
  mu_hybrid <- (effect_factor + 1)*((mu_female + mu_male)/2)
  
  # --- Norm Factors (Updated for Replicates) ---
  # Generate a matrix of norm factors: [n_families x n_reps]
  norm_generator <- function(nf, n_f, n_r) {
    if(length(nf) == 2) return(matrix(runif(n_f * n_r, min=nf[1], max=nf[2]), nrow=n_f))
    return(matrix(1, nrow=n_f, ncol=n_r))
  }
  F_NF_Mat <- norm_generator(female_norm_factor, n_families, n_reps)
  M_NF_Mat <- norm_generator(male_norm_factor, n_families, n_reps)
  H_NF_Mat <- norm_generator(hybrid_norm_factor, n_families, n_reps)
  
  # --- Dispersion Generation (Updated for Scenarios) ---
  mp <- (mu_female + mu_male)/2
  mp_seeds <- apply(mp, 1, median)
  
  if(MP_DISP_func == "low"){ 
    mp_seeding_disps <- sapply(mp_seeds, inverse_func)
  } else if(MP_DISP_func == "high"){
    # Assumes sim_params_med is loaded in global env
    mp_seeding_disps <- get_sim_dispersion(mp_seeds, sim_params_med, scale_noise=1)
  } else if(MP_DISP_func == "extreme") {
    mp_seeding_disps <- get_sim_dispersion(mp_seeds, sim_params_med, scale_noise=2)
  } else{
    stop("Dispersion is not valid!")
  }
  dispersion <- matrix(mp_seeding_disps, nrow = n_genes, ncol = n_families)
  dispersion[dispersion < 0.001] <- 0.001
  # print(summary(c(dispersion))) # Optional print
  
  # --- 2. STORAGE LISTS ---
  list_qval_2sLRT <- list() # Stores q-values for EVERY replicate of EVERY run
  list_qval_Rep   <- list() # Stores q-values for EVERY run (using all reps)
  
  # Helper to generate 3D Lambda
  calc_lambda_3D <- function(mu_mat, nf_mat) {
    arr <- array(0, dim = c(nrow(mu_mat), ncol(mu_mat), ncol(nf_mat)))
    for(r in 1:ncol(nf_mat)) arr[,,r] <- sweep(mu_mat, 2, nf_mat[,r], `*`)
    return(arr)
  }
  Lambda_F <- calc_lambda_3D(mu_female, F_NF_Mat)
  Lambda_M <- calc_lambda_3D(mu_male, M_NF_Mat)
  Lambda_H <- calc_lambda_3D(mu_hybrid, H_NF_Mat)
  
  # Helper to sample counts
  gen_counts_3D <- function(lam_arr, disp_mat) {
    disp_vec <- rep(as.vector(disp_mat), dim(lam_arr)[3]) 
    vals <- rnegbin(length(lam_arr), mu=as.vector(lam_arr), theta=1/disp_vec)
    array(vals, dim=dim(lam_arr))
  }
  
  # --- 3. SIMULATION LOOP ---
  for(iter in 1:repetition) {
    if(iter %% 10 == 0) cat(sprintf("Run %d/%d...\n", iter, repetition))
    
    # A. Generate Data
    Arr_F <- gen_counts_3D(Lambda_F, dispersion)
    Arr_M <- gen_counts_3D(Lambda_M, dispersion)
    Arr_H <- gen_counts_3D(Lambda_H, dispersion)
    
    # ----------------------------------------------------------------------
    # METHOD A: 2sLRT (Applied to EACH replicate individually)
    # ----------------------------------------------------------------------
    for(r in 1:n_reps) {
      # Extract Single Replicate Slice
      mat_F <- Arr_F[,,r]; mat_M <- Arr_M[,,r]; mat_H <- Arr_H[,,r]
      
      rownames(mat_F) <- rownames(mat_M) <- rownames(mat_H) <- paste0("Gene", 1:n_genes)
      colnames(mat_F) <- colnames(mat_M) <- colnames(mat_H) <- paste0("Family", 1:n_families)
      
      # Run Standard 2sLRT (Normalization = TRUE to estimate factors from slice)
      res_unrep <- twoStageLRT(hybrid_mat=mat_H, female_mat=mat_F, male_mat=mat_M,
                               disp_mat=NULL, normalization=TRUE, norm_factors=NULL, # use DESeq normalization and dispersion estimation
                               n_fam_thres=n_fam_thres, gene_groups_no=gene_groups_no)
      
      # Accumulate results
      list_qval_2sLRT[[length(list_qval_2sLRT)+1]] <- res_unrep$qval_NB
    }
    cat("For 2sLRT analysis: iteration through all reps done \n")
    
    # ----------------------------------------------------------------------
    # METHOD B: REPLICATED (Using ALL replicates)
    # ----------------------------------------------------------------------
    
    # 1. Flatten Data [Genes x (Families * Reps)]
    # Order: Fam1(R1..Rn), Fam2(R1..Rn)...
    mat_F_flat <- matrix(Arr_F, nrow=n_genes) 
    mat_M_flat <- matrix(Arr_M, nrow=n_genes)
    mat_H_flat <- matrix(Arr_H, nrow=n_genes)
    Big_Counts <- cbind(mat_F_flat, mat_M_flat, mat_H_flat)
    
    # 2. Global Dispersion Estimation
    # Group: Must define Family+Variety combo
    fam_indices <- rep(1:n_families, times=n_reps)
    Big_Group <- factor(c(paste0("F", fam_indices, "F"), 
                          paste0("F", fam_indices, "M"), 
                          paste0("F", fam_indices, "H")))
    
    # Estimate
    Big_C <- estimateSizeFactorsForMatrix(Big_Counts); if(any(is.na(Big_C))) Big_C[] <- 1
    y_dge <- DGEList(counts=Big_Counts, group=Big_Group)
    cat("For replicated analysis: size factor estimation done \n")
    ##### use all samples to estimate dispersion is too time consuming
    # y_dge$offset <- matrix(log(Big_C), nrow=nrow(Big_Counts), ncol=ncol(Big_Counts), byrow=TRUE)
    # y_dge <- estimateDisp(y_dge, model.matrix(~0+Big_Group), robust=TRUE)
    # global_phi <- y_dge$tagwise.dispersion
    ##### use the first 20 families for dispersion estimaiton
    n_sub <- min(20, n_families) # Use 20 families (or all if <20)
    # The matrix columns are: [F block] [M block] [H block]
    # fam_indices corresponds to ONE block. We repeat the mask 3 times.
    mask_block <- fam_indices <= n_sub
    keep_cols <- c(mask_block, mask_block, mask_block)
    y_sub <- DGEList(counts=Big_Counts[, keep_cols], group=droplevels(Big_Group[keep_cols]))
    y_sub$offset <- matrix(log(Big_C[keep_cols]), nrow=n_genes, ncol=sum(keep_cols), byrow=TRUE)
    y_sub <- estimateDisp(y_sub, model.matrix(~0+y_sub$samples$group), robust=TRUE)
    global_phi <- y_sub$tagwise.dispersion
    cat("For replicated analysis: dispersion estimation done \n")
    
    # 3. Family-wise Testing
    p_mat_rep <- matrix(NA, nrow=n_genes, ncol=n_families)
    
    # Helpers for indexing
    n_cols <- n_families * n_reps
    C_F_vec <- Big_C[1:n_cols]; C_M_vec <- Big_C[(n_cols+1):(2*n_cols)]; C_H_vec <- Big_C[(2*n_cols+1):(3*n_cols)]
    
    for(fam in 1:n_families) {
      cols <- seq(from = fam, by = n_families, length.out = n_reps)
      
      # Extract data for this family
      y_f <- mat_F_flat[, cols, drop=FALSE]
      y_m <- mat_M_flat[, cols, drop=FALSE]
      y_h <- mat_H_flat[, cols, drop=FALSE]
      c_f <- C_F_vec[cols]; c_m <- C_M_vec[cols]; c_h <- C_H_vec[cols]
      
      # Run Newton-Raphson LRT gene-by-gene
      fam_pvals <- numeric(n_genes)
      for(g in 1:n_genes) {
        # Using your Newton-Raphson helper defined outside
        fam_pvals[g] <- test_Family_LRT(list(y_f[g,], y_m[g,], y_h[g,]), list(c_f, c_m, c_h), global_phi[g])
      }
      p_mat_rep[, fam] <- fam_pvals
    }
    
    # Accumulate results
    list_qval_Rep[[iter]] <- apply(p_mat_rep, 2, function(x){p.adjust(x, method = "fdr")})
    cat("For replicated analysis: testing done \n")
    
  }
  
  # ============================================================================
  # 4. EVALUATION & PLOTTING (Matched to your original style)
  # ============================================================================
  
  if(mean(effect_factor == 0) != 1){
    
    # A. FDR and TPR Curves (using your alpha_func)
    # alpha_func automatically averages over the list length
    alpha_res_2sLRT <- alpha_func(list_qval_2sLRT, effect_factor) 
    alpha_res_Rep   <- alpha_func(list_qval_Rep, effect_factor)
    alphas <- alpha_res_2sLRT$alphas
    
    # FDR Plot Data
    FDR_avg_2sLRT <- alpha_res_2sLRT$FDR_avg
    FDR_avg_Rep   <- alpha_res_Rep$FDR_avg
    
    cat(paste("\n FDR at alpha=0.05 is, 2sLRT=", round(FDR_avg_2sLRT[alphas == 0.05],4),
              ", Rep=", round(FDR_avg_Rep[alphas == 0.05],4), sep =""))
    
    # TPR Plot Data
    TPR_avg_2sLRT <- alpha_res_2sLRT$TPR_avg
    TPR_avg_Rep   <- alpha_res_Rep$TPR_avg
    
    cat(paste("\n TPR at alpha=0.05 is, 2sLRT=", round(TPR_avg_2sLRT[alphas == 0.05],4),
              ", Rep=", round(TPR_avg_Rep[alphas == 0.05],4), sep =""))
    
    # B. ROC Curves (using your mean_ROC_func)
    roc_2sLRT <- mean_ROC_func(list_qval_2sLRT, effect_factor, direction = ">")
    roc_Rep   <- mean_ROC_func(list_qval_Rep, effect_factor, direction = ">")
    
    cat(paste("\n 2sLRT: pAUC =", roc_2sLRT$partial_AUC_avg, " SE =", roc_2sLRT$partial_AUC_SE, "\n"))
    cat(paste("Rep: pAUC =", roc_Rep$partial_AUC_avg, " SE =", roc_Rep$partial_AUC_SE, "\n"))
    
    # --- PLOTTING ---
    
    # FDR Plot
    p_FDR_avg <- ggplot(data = data.frame(alpha = rep(alphas,2), 
                                          FDR = c(FDR_avg_2sLRT, FDR_avg_Rep), 
                                          Method = rep(c("2sLRT","Replicated"), each = length(alphas))))+
      geom_line(aes(x = alpha, y = FDR, color = Method), linewidth = 1)+
      geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
      xlim(0, 0.2)+
      labs(title = "", y="Average FDR", x="Nominal FDR", color = "")+
      theme_bw()+
      theme(legend.position = "bottom")
    
    # TPR Plot
    p_TPR_avg <- ggplot(data = data.frame(alpha = rep(alphas,2), 
                                          TPR = c(TPR_avg_2sLRT, TPR_avg_Rep), 
                                          Method = rep(c("2sLRT","Replicated"), each = length(alphas))))+
      geom_line(aes(x = alpha, y = TPR, color = Method), linewidth = 1)+
      xlim(0, 0.2)+
      ylim(0, 1)+
      labs(title = "", y="Average TPR", x="Nominal FDR", color = "")+
      theme_bw()+
      theme(legend.position = "bottom")
    
    # ROC Plot
    p_roc_avg <- ggplot(data = data.frame(TPR = c(roc_2sLRT$tpr_actual, roc_Rep$tpr_actual), 
                                          FPR = rep(roc_2sLRT$fpr_actual, 2), 
                                          Method = factor(rep(c("2sLRT","Replicated"), each = length(roc_2sLRT$fpr_actual)),
                                                          levels=c("2sLRT","Replicated"))))+
      geom_line(aes(x = FPR, y = TPR, color = Method, linetype = Method), linewidth = 1)+
      xlim(0, 0.2)+
      ylim(0, 1)+
      labs(title = "", color="", linetype = "")+
      theme_bw()+
      theme(legend.position = "bottom")
    
    # Combine
    finalp <- ggarrange(p_FDR_avg, p_TPR_avg, p_roc_avg, nrow = 1, common.legend = TRUE, legend = "bottom")
    
  } else {
    finalp <- NULL
  }
  
  if(!is.null(output_file)){
    ggsave(paste(output_file, ".pdf", sep = ""), finalp, width = 12, height = 4, units = "in")
  }
  
  return(finalp)
}

# ==============================================================================
# 6. EXECUTION BLOCK (Argument Parsing)
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0) {
  pattern_arg <- args[1]
  N0_arg <- as.numeric(args[2])
  K_arg <- as.numeric(args[3])
  ngenes_arg <- as.numeric(args[4])
  nfamilies_arg <- as.numeric(args[5])
  disp_arg <- args[6]
  n_reps_arg <- as.numeric(args[7]) 
  output_filename_arg <- args[8]
  
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
  
  # case <- list(n_genes = 100, n_families = 100,
  #              sig_pos = 3, sig_neg = -3/4,
  #              rate_genes = rate_genes,
  #              mu_female = mu_female,
  #              mu_male = mu_male,
  #              effect_factor = effect_factor,
  #              female_norm_factor = c(0.8, 1.3), male_norm_factor = c(0.8, 1.3), hybrid_norm_factor = c(0.8, 1.3),
  #              MP_DISP_func = "low",
  #              n_fam_thres =5, gene_groups_no = 1, repetition = 3, output_file = "aaa", n_reps = 3)

  case <- list(n_genes = ngenes_arg, n_families = nfamilies_arg,
               sig_pos = 3, sig_neg = -3/4,
               rate_genes = rate_genes,
               mu_female = mu_female,
               mu_male = mu_male, 
               effect_factor = effect_factor,
               female_norm_factor = c(0.8, 1.3), male_norm_factor = c(0.8, 1.3), hybrid_norm_factor = c(0.8, 1.3),
               MP_DISP_func = disp_arg, 
               n_fam_thres = N0_arg, gene_groups_no = K_arg, repetition = 10, 
               n_reps = n_reps_arg,
               output_file = output_filename_arg)
  
  print(case)
  if(is.character(case$rate_genes)){case$rate_genes <- eval(parse(text = case$rate_genes))}
  if(is.character(case$mu_female)){case$mu_female <- eval(parse(text = case$mu_female))}
  if(is.character(case$mu_male)){case$mu_male <- eval(parse(text = case$mu_male))}
  if(is.character(case$effect_factor)){case$effect_factor <- eval(parse(text = case$effect_factor))}
   res <- Heter_pip_Comparison(n_genes = case$n_genes, n_families = case$n_families,
                                 mu_female = case$mu_female, mu_male = case$mu_male,
                                 effect_factor = case$effect_factor,
                                 female_norm_factor = case$female_norm_factor,
                                 male_norm_factor = case$male_norm_factor,
                                 hybrid_norm_factor = case$hybrid_norm_factor,
                                 MP_DISP_func = case$MP_DISP_func,
                                 n_fam_thres = case$n_fam_thres,
                                 gene_groups_no = case$gene_groups_no,
                                 repetition = case$repetition, 
                                 n_reps = case$n_reps,
                                 output_file = case$output_file)
}
