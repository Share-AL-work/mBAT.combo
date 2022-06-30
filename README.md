
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mBAT

<!-- badges: start -->
<!-- badges: end -->

Multivariate set-Based Association Test (mBAT), a novel gene-based
method, to improve power of identifying genes harbouring variants with
masking effects. We show by extensive simulations that mBAT has greater
power than the sum-χ^2 strategy (taking fastBAT and MAGMA as examples)
in the presence of masking effects. To maximise overall power regardless
of the masking effects, we further proposed a hybrid method, mBAT-combo,
by combining mBAT and fastBAT test statistics through a Cauchy
combination method, a recently developed method to combine different
test statistics without knowing the correlation structure(Liu Y, Chen S,
Li Z, Morrison AC, Boerwinkle E, Lin X. ACAT: A Fast and Powerful p
Value Combination Method for Rare-Variant Analysis in Sequencing
Studies. Am J Hum Genet. 2019 Mar 7;104(3):410-421. doi:
10.1016/j.ajhg.2019.01.002. PMID: 30849328; PMCID: PMC6407498.).

## Installation

You can install the development version of mBAT from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Share-AL-work/mBAT.combo")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(mBAT)
```

## 1. set path

``` r
path <- system.file("extdata",package = "mBAT")
list.files(path)
#> [1] "1000G_eur_unrel_ukbcom22.bim" "Chr22.fastGWA.ma"            
#> [3] "fb_Chr22_1kg.gene.fastbat"    "hg19_v40_glist_Pat_v1.txt"   
#> [5] "LD_Chr22"                     "mBAT_R_Chr22.txt"

bim_file <- paste0(path,"/","1000G_eur_unrel_ukbcom22.bim")
map_file <- paste0(path,"/","hg19_v40_glist_Pat_v1.txt")
assoc_file <- paste0(path,"/","Chr22.fastGWA.ma")
LD_path <- paste0(path,"/","LD_Chr22")
#result_path <- paste0(path,"/","extdata")
fastBAT_output <- paste0(path,"/","fb_Chr22_1kg.gene.fastbat")
```

## 2. Do analysis

``` r
prop=0.9
#res <- mBAT::mBAT_combo(bim_file,
#                  map_file,
#                  assoc_file,
#                  LD_path,
#                  result_path,
#                  fastBAT_output,
#                  prop,
#                  gene_annotate)
#head(res)
```

\#n that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
