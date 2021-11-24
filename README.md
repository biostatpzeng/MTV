# MTV: Mixed TWAS and mediated Variance estimation
========================================================================================================
## Background
Although several integrative methods have been proposed, how to incorporate eQTL more efficiently remains less understood. The mixed effects score test (MiST) was developed for SNP-set association studies by modeling impacts of genetic variants based on functional annotations while allowing for variant-specific influence. However, MiST conducts score-based test, thus cannot quantify the relative contribution of eQTL and genetic loci to phenotypic variation. On the other hand, as gene expressions are often directly unavailable or unmeasured in GWAS due to cost and unavailability of specimen, to investigate the relationship between unmeasured gene expression and trait/disease, transcriptome-wide association studies (TWASs) were developed to bridge such gap by imputing/predicting unmeasured gene expression via an external transcriptome reference panel, where both expressions and genotypes were available. Existing TWAS methods (e.g., prediXcan) make a relatively strong modeling assumption that cis-SNPs of a gene do not exhibit any direct effects (i.e., horizontal pleiotropy), in contrast to empirical evidence that genetic variants can influence phenotype both directly and indirectly. TWAS can be viewed as a special Mendelian randomization (MR) method, which assumes the horizontal pleiotropy is absent for used instruments (i.e., direct cis-SNP effects are zero); the failure of handling cis-SNP direct effects properly could lead to spurious associations in TWAS where some observed association signs do not necessarily represent the true relationship between genes and phenotype. Instead, such associations might be driven by cis-SNPs alone, which also partly explains why many identified genes in previous TWASs were generally located near or within clusters of associated GWAS loci.

To integrate eQTL mapping study into GWAS, we here proposed a novel statistical method, called **MTV (Mixed TWAS and mediated Variance estimation)**, by modeling the effects of cis-SNPs of a given gene as a function of eQTL. As would be shown, MTV formulates the eQTL integrative method (e.g., MiST) and individual-level TWAS methods (e.g., prediXcan) within a unified framework through mixed models and includes many prior methods/tests as special cases. We further justified MTV from another two statistical perspectives including mediation analysis and two-stage MR. To efficiently estimate unknown parameters in MTV, we developed a parameter expansion expectation maximum (PX-EM) algorithm that can be scalable to large-scale biobank data. We defined two useful quantities to measure relative genetic contributions of gene expression and its cis-SNPs to phenotypic variance. In addition, a fast and powerful likelihood ratio test method was proposed in MTV to jointly test the total effects of cis-SNPs and gene expression on the phenotype. With extensive simulation studies, we demonstrated that MTV can correctly maintain the type I error rate control when jointly testing the total genetic effects and was often more powerful to identify true association signals across various scenarios compared to existing methods. 

**[MTV](https://github.com/biostatpzeng/HMAT/blob/main/HMAT_function.R)** is implemented in R statistical environment.

## Example
For GWAS with individual genotyps and phenotype
```ruby
source("HMAT_function.R")
y <- read.table("y.txt",sep=""),head=F)[,1]
G2 <- read.table("snp_gwas.txt",head=F)
weight <- matrix(runif(m*7),m,7)

# Here, we assume, for simplicity, that these simulated weights are estimated from seven various gene expression
# prediction models. Then, actually, there are seven various TWAS analyses. For each TWAS, we can obtain its p value
# to evaluate the significance of the gene. Finally, we combine these p values into a single one using HMAT.

HMAT_individual(y,G2,weight,outcome="B")

$p_HMAT
[1] 0.9125274

$p_TWAS
[1] 0.9185220 0.8028170 0.7293027 0.9295604 0.7424362 0.9007160 0.9363431
```

## Cite
Ting Wang, Jiahao Qiao, Shuo Zhang, Yongyue Wei and [Ping Zeng](https://github.com/biostatpzeng) (2021). Simultaneous test and estimation of total genetic effect in eQTL integrative analysis through mixed models.

## Contact
We are very grateful to any questions, comments, or bugs reports; and please contact [Ping Zeng](https://github.com/biostatpzeng) via zpstat@xzhmu.edu.cn.

## Update
2021-11-24 MTV version 1.0.
