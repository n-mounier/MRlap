
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

:information\_source: `MRlap` is still under active development.

## Overview

`MRlap` is an R-package to perform two-sample Mendelian Randomisation
(MR) analyses for (potentially) overlapping samples, using only GWAS
summary statistics. MR estimates can be subject to different types of
biases due to the overlap between the exposure and outcome samples, the
use of weak instruments and Winner’s curse. Our approach simultaneously
accounts and corrects for all these biases. Estimating the corrected
effect using our approach can be performed as a sensitivity analysis: if
the corrected effect do not significantly differ from the observed
effect, then IVW-MR estimate can be safely used. However, when there is
a significant difference, corrected effects should be preferred as they
should be less biased, independently of the sample overlap.  
This package builds up on the
[`GenomicSEM`](https://github.com/GenomicSEM/GenomicSEM/) R-package to
perform cross-trait LD-score regression (LDSC) and the [`TwoSampleMR`]()
R-package for inverse-variance weighted (IVW-)MR analysis (and
instruments pruning).

There is only one function available:

  - **`MRlap()`**  
    main function that performs LDSC, IVW-MR analysis and provides a
    corrected causal effect estimate.

<!--More details about their usage can be found in the [manual](doc/bGWAS-manual.pdf).-->

## Installation

You can install the current version of `MRlap` with:

``` r
# Directly install the package from github
# install.packages("remotes")
##remotes::install_github("n-mounier/MRlap")
##library(MRlap)
```

<!--- Note: using remotes instead of devtools leads to re-build the package
and apparently, it may be a problem with R 3.4 and macOS, 
see https://stackoverflow.com/questions/43595457/alternate-compiler-for-installing-r-packages-clang-error-unsupported-option/43943631#43943631 --->

## Usage

To run the analysis with `MRlap` different inputs are needed: \#\#\#\#
1. The exposure and outcome GWAS summary statistics (`exposure` &
`outcome`): Can be a regular (space/tab/comma-separated) file or a
gzipped file (.gz) or a `data.frame`. Must contain the following
columns, which can have alternative names (not case sensitive):  

<ul>

SNP-identifier: `rs` or `rsid`, `snp`, `snpid`, `rnpid`  
Alternate (effect) allele: `a1` or `alt`, `alts`  
Reference allele: `a2` or `a0`, `ref`  
Z-statistics: `Z` or `zscore`  
Sample size: `N`

</ul>

If the Z-statistics is not present, it will be automatically calculated
from effect size and standard error, in which case the following columns
should be provided:  

<ul>

Effect-size: `b` or `beta`, `beta1`  
Standard error: `se` or `std`

</ul>

#### 2\. The input files for LDSC (`ld` & `hm3`):

These are needed by the
[`GenomicSEM`](https://github.com/GenomicSEM/GenomicSEM/) R-package.  
\- ld: \> Expects LD scores formated as required by the original LD
score regression software. Weights for the european population can be
obtained by downloading the eur\_w\_ld\_chr folder in the link below
(Note that these are the same weights provided by the original
developers of LDSC):
<https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v> - hm3 \> We
suggest using an (UNZIPPED) file of HAPMAP3 SNPs with some basic
cleaning applied (e.g., MHC region removed) that is supplied and created
by the original LD score regression developers and available here:
<https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2>

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
