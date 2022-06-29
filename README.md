
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
#library(mBAT)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
