% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Pois_3mean_LRT.R
\name{Pois_3mean_LRT}
\alias{Pois_3mean_LRT}
\title{Poisson Likelihood Ratio Test for Equal Family Mean
This function conducts Poisson LRT test for detection of null-family in a single family and for a single gene with single replication per genotype.
Denote the mean for gene g and female parent as mu_g1, male parent as mu_g2 and hybrid as mu_g3, the null hypothesis is Hg0: mu_g1 = mu_g2 = mu_g3.}
\usage{
Pois_3mean_LRT(y_g1, y_g2, y_g3, C_1 = 1, C_2 = 1, C_3 = 1)
}
\arguments{
\item{y_g1}{raw count for female.}

\item{y_g2}{raw count for male.}

\item{y_g3}{raw count for hybrid.}

\item{C_1}{normalization factor for female sample.}

\item{C_2}{normalization factor for male sample.}

\item{C_3}{normalization factor for hybrid sample.}
}
\value{
p-value.
}
\description{
Poisson Likelihood Ratio Test for Equal Family Mean
This function conducts Poisson LRT test for detection of null-family in a single family and for a single gene with single replication per genotype.
Denote the mean for gene g and female parent as mu_g1, male parent as mu_g2 and hybrid as mu_g3, the null hypothesis is Hg0: mu_g1 = mu_g2 = mu_g3.
}
\examples{
library(TwoStageLRT)
Pois_3mean_LRT(y_g1 = 10, y_g2 = 15, y_g3 = 20, C_1 = 1, C_2 = 1.2, C_3 = 0.8)
}
