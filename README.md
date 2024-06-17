
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TwoStageLRT

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

The goal of TwoStageLRT (2sLRT) is to detect gene expression mid-parent
heterosis (MPH) from unreplicated multi-family RNA-seq experiment.

If you want to use 2sLRT, you might have gene expression data for
multiple genes and multiple families, each family has three genotypes,
female parent, male parent and their hybrid genotype. The key feature of
using 2sLRT is there is only one sample for each genotype, i.e.,
unreplicated.

This test has two stages to deal with the problem of lacking
replications. In the first stage, dispersion paremters are estimated. In
the second stage, likelihood ratio test based on negative binomial
distribution given the estimated dispersion parameters is conducted.

## Installation

You can install the development version of TwoStageLRT like so:

``` r
# install.packages("devtools")
devtools::install_github("yunhuiqistat/TwoStageLRT")
```

## Website

Check the vignette at <https://yunhuiqistat.github.io/TwoStageLRT/>.

## Usage

To use 2sLRT, you need to have at least three data matrices:

- female_mat: female expression count with genes in rows, samples in
  columns.

- male_mat: male expression count with genes in rows, samples in
  columns.

- hybrid_mat: hybrids expression count with genes in rows, samples in
  columns.

You can also specify the dispersion matrix, normalization factors if you
have prior knowledge.

- disp_mat if “Pois”, use Poisson LRT and only output Pois results; if
  NULL, use dispersion estimation output will have both NBLRT and
  PoisLRT results; else, use provided disp_mat for NB test, output will
  have both NBLRT and PoisLRT results.

- normalization logical, whether or not need normalization for the raw
  data matrices, if FALSE, treat the normalization factors all 1.

- norm_factors if normalization = TRUE and norm_factors = NULL, use
  DESeq normalization factors; if normalization = TRUE and norm_factors
  is a list contains three vectors female, male and hybrid normalization
  factors, use the provided norm_factors in tests.

Besides, you should also specify the hyper-parameters for dispersion
estimation:

- n_fam_thres genes with at least n_fam_thres null families will be used
  for dispersion estimation, otherwise, will use loess prediction to get
  dispersion, default is 20.

- gene_groups_no the number of clusters for genes, default is 3.

2sLRT will give results based on your choice:

- if disp_mat = “Pois”, return Poisson LRT for heterosis pvalue and
  adjusted pvalue from BH adjustment;

- if disp_mat = NULL or provided, return both NB and Poisson LRT for
  heterosis pvalue and qvalue and genewise dispersion, and the
  estimation of mu_gi, which is the female/male/hybrid_mat_normed.

``` r
library(TwoStageLRT)
data(sim_data)
res <- twoStageLRT(female_mat = sim_data$female_mat, male_mat = sim_data$male_mat,
                   hybrid_mat = sim_data$hybrid_mat, disp_mat = NULL,
                   normalization = TRUE, norm_factors = NULL,
                   n_fam_thres = 10, gene_groups_no = 2)
#> Use DESeq estimated normalization factors. 
#> Gene clusters size: 
#> jacc_clust
#>  1  2 
#> 38  6 
#> Dispersion estimation for cluster 1 the number of families is 10 , number of genes is  38  
#> Dispersion estimation for cluster 2 the number of families is 14 , number of genes is  6  
#> Dispersion estimation summary:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#> 0.02848 0.07723 0.09988 0.09786 0.11958 0.20517      56 
#> Output is from NB LRT for heterosis.
```

``` r
summary(unlist(c(res$pval_NB)))
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> 0.000000 0.006074 0.140645 0.287458 0.537483 0.999942
```

``` r
summary(unlist(c(res$pval_Pois)))
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> 0.000000 0.000000 0.000000 0.067466 0.002527 0.999684
```

``` r
summary(res$DISP)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.02848 0.08031 0.09892 0.09567 0.11058 0.20517
```
