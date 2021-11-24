# MTV: Mixed TWAS and mediated Variance estimation
==============================================================================================
## Background
Although several integrative methods have been proposed, how to incorporate eQTL more efficiently remains less understood. The mixed effects score test (MiST) was developed for SNP-set association studies by modeling impacts of genetic variants based on functional annotations while allowing for variant-specific influence. However, MiST conducts score-based test, thus cannot quantify the relative contribution of eQTL and genetic loci to phenotypic variation. On the other hand, as gene expressions are often directly unavailable or unmeasured in GWAS due to cost and unavailability of specimen, to investigate the relationship between unmeasured gene expression and trait/disease, transcriptome-wide association studies (TWASs) were developed to bridge such gap by imputing/predicting unmeasured gene expression via an external transcriptome reference panel, where both expressions and genotypes were available. Existing TWAS methods (e.g., prediXcan) make a relatively strong modeling assumption that cis-SNPs of a gene do not exhibit any direct effects (i.e., horizontal pleiotropy), in contrast to empirical evidence that genetic variants can influence phenotype both directly and indirectly. TWAS can be viewed as a special Mendelian randomization (MR) method, which assumes the horizontal pleiotropy is absent for used instruments (i.e., direct cis-SNP effects are zero); the failure of handling cis-SNP direct effects properly could lead to spurious associations in TWAS where some observed association signs do not necessarily represent the true relationship between genes and phenotype. Instead, such associations might be driven by cis-SNPs alone, which also partly explains why many identified genes in previous TWASs were generally located near or within clusters of associated GWAS loci.

To integrate eQTL mapping study into GWAS, we here proposed a novel statistical method, called **MTV (Mixed TWAS and mediated Variance estimation)**, by modeling the effects of cis-SNPs of a given gene as a function of eQTL. As would be shown, MTV formulates the eQTL integrative method (e.g., MiST) and individual-level TWAS methods (e.g., prediXcan) within a unified framework through mixed models and includes many prior methods/tests as special cases. We further justified MTV from another two statistical perspectives including mediation analysis and two-stage MR. To efficiently estimate unknown parameters in MTV, we developed a parameter expansion expectation maximum (PX-EM) algorithm that can be scalable to large-scale biobank data. We defined two useful quantities to measure relative genetic contributions of gene expression and its cis-SNPs to phenotypic variance. In addition, a fast and powerful likelihood ratio test method was proposed in MTV to jointly test the total effects of cis-SNPs and gene expression on the phenotype. With extensive simulation studies, we demonstrated that MTV can correctly maintain the type I error rate control when jointly testing the total genetic effects and was often more powerful to identify true association signals across various scenarios compared to existing methods. 

**[MTV](https://github.com/biostatpzeng/MTV)** is implemented in R statistical environment.

## Example
For the parameter estimation in MTV
```ruby
library(Rcpp)
library(RcppArmadillo)
sourceCpp("lmm_pxem.cpp")
sourceCpp("LRTsim.cpp")
source("LRTsim.R")

fit = lmm_pxem(y, X=cbind(1, E), G=snp, PXEM=TRUE, maxIter=1000)

sb = mean(snp*snp)
m = dim(snp)[2]
EB = E*fit$alpha[2]
sigmat2=c(fit$theta)[2]
sigmate=c(fit$theta)[1]
pve = (sb*sigmat2*m + var(EB))/(var(EB) + sb*sigmat2*m + sigmate)
pge =                 var(EB) /(var(EB) + sb*sigmat2*m)

//' @param y  response variable for GWAS data
//' @param X  covariates for GWAS data, here GReX should be included
//' @param G  normalized genotype (cis-SNPs) matrix for GWAS
//' @param maxIter  maximum iteration (default is 1000)

```
For joint effect test using LRT in MTV
```ruby
library(Rcpp)
library(RcppArmadillo)
sourceCpp("lmm_pxem.cpp")
sourceCpp("LRTsim.cpp")
source("LRTsim.R")

fit0 = lm(y~Z)
fit1 = lmm_pxem(y, X=cbind(1, Z, E), G=G, PXEM=TRUE, maxIter=1000)
simLike <- eLRT(Z = Z, E = E, G = G, nsim=1e5, parallel=c("multicore"), ncpus = 4L) ## exact LRT

obsLike = c((fit1$loglik - logLik(fit0))*2)
p1 = mean(simLike >= obsLike)
p2 = aLRT(simLike, c(1e3,1e4,1e5)) ## approximate LRT

//' @param y  response variable for GWAS data
//' @param Z  covariates for GWAS data
//' @param E  GReX
//' @param G  normalized genotype (cis-SNPs) matrix for GWAS

```

## Cite
Ting Wang, Jiahao Qiao, Shuo Zhang, Yongyue Wei and [Ping Zeng](https://github.com/biostatpzeng) (2021). Simultaneous test and estimation of total genetic effect in eQTL integrative analysis through mixed models.

## Contact
We are very grateful to any questions, comments, or bugs reports; and please contact [Ping Zeng](https://github.com/biostatpzeng) via zpstat@xzhmu.edu.cn.

## Update
2021-11-24 MTV version 1.0.
