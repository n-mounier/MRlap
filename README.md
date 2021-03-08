
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MRlap <img src="inst/Figures/logo.png" align="right" height=180/>

<!--- 
# https://github.com/GuangchuangYu/hexSticker
library(hexSticker)
imgurl <- "inst/Figures/MRlap.png"
sticker(imgurl, 
        package="", p_size=8, p_color="black",
        h_fill="black", h_color="black",
        s_x=0.95, s_y=1.05, s_width=0.95,
        filename="inst/Figures/logo.png", dpi=2000) --->

<!--- :arrow_right: ESHG poster is available [here]().  --->

:information\_source: `MRlap` is still under active development

## Overview

`MRlap` is an R-package …  
builds up on [`GenomicSEM`]() to perform cross-trait LD-score regression
and \[`TwoSampleMR`\] for IVW analysis (and IVs clumping).

There is only one function available:

  - **`MRlap()`**  
    main function that …

## Installation

You can install the current version of `MRlap` with:

``` r
# Directly install the package from github
# install.packages("remotes")
#remotes::install_github("n-mounier/MRlap")
#library(MRlap)
```

<!--- Note: using remotes instead of devtools leads to re-build the package
and apparently, it may be a problem with R 3.4 and macOS, 
see https://stackoverflow.com/questions/43595457/alternate-compiler-for-installing-r-packages-clang-error-unsupported-option/43943631#43943631 --->

## Usage

To run the analysis with `MRlap` different inputs are needed:

### Results

##### Aditionnaly, if `save_logfiles=TRUE`, LDSC log files are created in the current working directory :

  - **<exposure>.log** - log file

## Runtime

<!---Analysis using all the 38 prior GWASs available, for a conventional GWAS containing ~7M SNPs in common with the prior studies ~ 25 minutes.

Analysis using 6 prior GWASs, for a conventional GWAS containing ~ 300,000 SNPs in common with prior studies (see example A) ~ 2 minutes.--->

<font color="grey"><small> Results from analyses performed on a MacBook
Pro (Early 2015) - Processor : 2.9 GHz Intel Core i5 - Memory : 8 GB
1867 MHz DDR3.</font> </small>

## Citation

If you use the `MRlap` package, please cite:

<!--- [Ninon Mounier, Zoltán Kutalik, bGWAS: an R package to perform Bayesian Genome Wide Association Studies, Bioinformatics](https://doi.org/10.1093/bioinformatics/btaa549) --->

## Contact

<mounier.ninon@gmail.com>
