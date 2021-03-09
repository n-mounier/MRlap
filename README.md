
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MRlap <img src="inst/Figures/logo.png" align="right" height=180/>

<!--- 
# https://github.com/GuangchuangYu/hexSticker
library(hexSticker)
imgurl <- "inst/Figures/MRlap2.png"
sticker(imgurl, 
        package="", p_size=8, p_color="black",
        h_fill="black", h_color="black",
        s_x=0.95, s_y=0.90, s_width=0.95, # MRlap2
        # s_x=0.95, s_y=1.05, s_width=0.95, # MRlap
        filename="inst/Figures/logo.png", dpi=2000) --->

<!--- :arrow_right: ESHG/EMGM?? poster is available [here]().  --->

:information\_source: `MRlap` is still under active development.

## Overview

`MRlap` is an R-package to perform two-sample Mendelian Randomisation
(MR) analyses using (potentially) overlapping samples, relying only on
GWAS summary statistics. MR estimates can be subject to different types
of biases due to the overlap between the exposure and outcome samples,
the use of weak instruments and Winner’s curse. Our approach
simultaneously accounts and corrects for all these biases, using
cross-trait LD-score regression (LDSC) to approximate the overlap.
Estimating the corrected effect using our approach can be performed as a
sensitivity analysis: if the corrected effect do not significantly
differ from the observed effect, then IVW-MR estimate can be safely
used. However, when there is a significant difference, corrected effects
should be preferred as they should be less biased, independently of the
sample overlap.  
This package builds up on the
[`GenomicSEM`](https://github.com/GenomicSEM/GenomicSEM/) R-package to
perform cross-trait LDSC and the
[`TwoSampleMR`](https://github.com/MRCIEU/TwoSampleMR/) R-package for
inverse-variance weighted (IVW-)MR analysis (and instruments pruning).

There is only one function available:

  - **`MRlap()`**  
    main function that performs LDSC, IVW-MR analysis and provides a
    corrected causal effect estimate.

More details about its usage can be found in the
[manual](doc/MRlap-manual.pdf).

## Installation

You can install the current version of `MRlap` with:

``` r
# Directly install the package from github
# install.packages("remotes")
remotes::install_github("n-mounier/MRlap")
library(MRlap)
```

<!--- Note: using remotes instead of devtools leads to re-build the package
and apparently, it may be a problem with R 3.4 and macOS, 
see https://stackoverflow.com/questions/43595457/alternate-compiler-for-installing-r-packages-clang-error-unsupported-option/43943631#43943631 --->

## Usage

To run the analysis with `MRlap` different inputs are needed:

#### 1\. The exposure and outcome GWAS summary statistics (`exposure` & `outcome`):

Can be a regular (space/tab/comma-separated) file or a gzipped file
(.gz) or a `data.frame`. Must contain the following columns, which can
have alternative names (not case sensitive):  

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

  - ld:

> Expects LD scores formated as required by the original LD score
> regression software. Weights for the european population can be
> obtained by downloading the eur\_w\_ld\_chr folder in the link below
> (Note that these are the same weights provided by the original
> developers of LDSC):
> <https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v>

  - hm3:

> We suggest using an (UNZIPPED) file of HAPMAP3 SNPs with some basic
> cleaning applied (e.g., MHC region removed) that is supplied and
> created by the original LD score regression developers and available
> here:
> <https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2>

### Analysis

Before running the examples, please make sure to have downloaded the
input files for LDSC. You may also need to modify the `ld` & `hm3`
parameters to indicate the correct paths.

  - **Example A**

<!-- end list -->

``` r
# Using ~100K samples for BMI/SBP, with 0% of sample overlap
# Note that here the overlap is known (since we generated the data) but the MRlap
# function works even the overlap is unkown (overlap is *not* a parameter of the function) 
# as it uses cross-trait LDSC to approximate it
# (only weak-instrument bias and Winner's curse)
# (1,150,000 SNPs - stored in gzipped files)
BMI <- system.file("data/", "BMI_Data.tsv.gz", package="MRlap")
SBP <- system.file("data/", "SBP_Data.tsv.gz", package="MRlap")

# MR instruments will be selected using default parameter (5e-8) and distance-pruned (500Kb),
# No file will be saved.
A = MRlap(exposure = SmallExposure_Data,
          exposure_name = "simulated_exposure",
          outcome = SmallOutcome_Data,
          outcome_name = "simulated_outcome",
          ld = "~/eur_w_ld_chr",
          hm3 = "~/w_hm3.noMHC.snplist")
```

<details>

<summary>Show log</summary>

    ```
    ## <<< Preparation of analysis >>> 
    ##  > Checking parameters 
    ##  The p-value threshold used for selecting MR instruments is: 5e-08 
    ##  The distance used for pruning MR instruments is:  500 Kb  
    ##  > Processing exposure (BMI_100Ksample) summary statistics...  
    ##  # Preparation of the data...  
    ##  The data.frame used as input is: "BMI_Data.tsv.gz".   
    ##     SNPID column, ok - CHR column, ok - POS column, ok - ALT column, ok - REF column, ok - Z column, ok - N column, ok  
    ##  > Processing outcome (SBP_100Ksample) summary statistics...  
    ##  # Preparation of the data...  
    ##  The data.frame used as input is: "SBP_Data.tsv.gz".   
    ##     SNPID column, ok - CHR column, ok - POS column, ok - ALT column, ok - REF column, ok - Z column, ok - N column, ok  
    ##  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
    ##  <<< Performing cross-trait LDSC >>>   
    ##  > Munging exposure data...  
    ##  > Munging outcome data...  
    ##  > Running cross-trait LDSC...  
    ##    Please consider saving the log files and checking them to ensure that all columns were interpreted correctly and no warnings were issued for any of the summary statistics files 
    ##  > Cleaning temporary files...  
    ##  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
    ##  <<< Running IVW-MR >>>   
    ##  > Identifying IVs...  
    ##      668 IVs with p < 5e-08  
    ##      0 IVs excluded - more strongly associated with the outcome than with the exposure, p < 1e-03  
    ##     Pruning : distance :  500 Kb  
    ##      39 IVs left after pruning  
    ##  > Performing MR  
    ##      IVW-MR observed effect: 0.0856 ( 0.0398 ) 
    ##  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
    ##  <<< Estimating corrected effect >>>   
    ##  > Estimating genetic architecture parameters...  
    ##  > Estimating corrected effect...  
    ##      corrected effect: 0.115 ( 0.0535 ) 
    ##      covariance between observed and corrected effect: 0.00215   
    ##  > Testing difference between observed and corrected effect...  
    ##  Runtime of the analysis:  4  minute(s) and  36  second(s).
    ```

</details>

  - **Example B**

<!-- end list -->

``` r
# Using simulated exposure/outcome data 
# standard settings scenario, with 100% of sample overlap
# Note that here the overlap is known (since we generated the data) but the MRlap
# function works even the overlap is unkown (overlap is *not* a parameter of the function) 
# as it uses cross-trait LDSC to approximate it
# (~400,000 SNPs - stored as data.frames)

# MR instruments will be selected using a more stringent threshold (5e-10) and LD-pruned (500Kb - r2=0.05),
# No file will be saved.
B = MRlap(exposure = SmallExposure_Data,
          exposure_name = "simulated_exposure",
          outcome = SmallOutcome_Data,
          outcome_name = "simulated_outcome",
          ld = "~/eur_w_ld_chr",
          hm3 = "~/w_hm3.noMHC.snplist",
          MR_threshold = 5e-10,
          MR_pruning_LD = 0.05)
```

<details>

<summary>Show log</summary>

    ```
    ## <<< Preparation of analysis >>> 
    ##  > Checking parameters 
    ##  The p-value threshold used for selecting MR instruments is: 5e-10 
    ##  The distance used for pruning MR instruments is:  500 Kb  
    ##  The LD threshold used for pruning MR instruments is: 0.05 
    ##  > Processing exposure (simulated_exposure) summary statistics...  
    ##  # Preparation of the data...  
    ##  The data.frame used as input is: "SmallExposure_Data".   
    ##     SNPID column, ok - CHR column, ok - POS column, ok - ALT column, ok - REF column, ok - Z column, ok - N column, ok  
    ##  > Processing outcome (simulated_outcome) summary statistics...  
    ##  # Preparation of the data...  
    ##  The data.frame used as input is: "SmallOutcome_Data".   
    ##     SNPID column, ok - CHR column, ok - POS column, ok - ALT column, ok - REF column, ok - Z column, ok - N column, ok  
    ##  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
    ##  <<< Performing cross-trait LDSC >>>   
    ##  > Munging exposure data...  
    ##  > Munging outcome data...  
    ##  > Running cross-trait LDSC...  
    ##    Please consider saving the log files and checking them to ensure that all columns were interpreted correctly and no warnings were issued for any of the summary statistics files 
    ##  > Cleaning temporary files...  
    ##  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
    ##  <<< Running IVW-MR >>>   
    ##  > Identifying IVs...  
    ##      106 IVs with p < 5e-10  
    ##      0 IVs excluded - more strongly associated with the outcome than with the exposure, p < 1e-03  
    ##     Pruning : distance :  500 Kb  - LD threshold :  0.05  
    ##      32 IVs left after pruning  
    ##  > Performing MR  
    ##      IVW-MR observed effect: 0.251 ( 0.0242 ) 
    ##  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
    ##  <<< Estimating corrected effect >>>   
    ##  > Estimating genetic architecture parameters...  
    ##  > Estimating corrected effect...  
    ##      corrected effect: 0.236 ( 0.0277 ) 
    ##      covariance between observed and corrected effect: 0.00067  
    ##  > Testing difference between observed and corrected effect...  
    ##  Runtime of the analysis:  3  minute(s) and  41  second(s).
    ```

</details>

### Results

**`MRlap()`** returns a named list containing the following results:

  - MRcorrection

<ul>

“observed\_effect” : IVW-MR observed causal effect estimate,  
“observed\_effect\_se” : IVW-MR observed causal effect estimate standard
error,  
“corrected\_effect” : corrected causal effect estimate,  
“corrected\_effect\_se” : corrected causal effect estimate standard
error,  
“test\_difference” : test statistic used to test for differences between
observed and corrected effects,  
“p\_difference” : p-value corresponding to the test statistic used to
test for differences between observed and corrected effects.

</ul>

  - LDSC

<ul>

“h2\_exp” : exposure heritability estimate,  
“h2\_exp\_se” : exposure heritability standard error,  
“int\_exp” : exposure intercept,  
“h2\_out” : outcome heritability estimate,  
“h2\_out\_se” : outcome heritability standard error,  
“int\_out” : outcome intercept,  
“gcov” : genetic covariance estimate,  
“gcov\_se” : genetic covariance estimate standard error,  
“rg” : genetic correlation estimate, “int\_crosstrait” : cross-trait
intercept estimate, “int\_crosstrait\_se”: cross-trait intercept
estimate standard error.

</ul>

  - GeneticArchitecture

<ul>

“polygenicity” : exposure polygenicity estimate,  
“perSNP\_heritability” : exposure per-SNP heritability estimate.

</ul>

##### Aditionnaly, if `save_logfiles=TRUE`, LDSC log files are created in the current working directory :

  - **<exposure_name>.log** - exposure cleaning/munging log file  

  - **<outcome_name>.log** - outcome cleaning/munging log file  

  - **<exposure_name>.sumstats.gz\_<outcome_name>.sumstats.gzldsc.log**
    - cross-trait LDSC log file

  - **Example A**

<!-- end list -->

``` r
# structure of the results
str(A)
```

    ## List of 3
    ##  $ MRcorrection       :List of 6
    ##   ..$ observed_effect    : num 0.0856
    ##   ..$ observed_effect_se : num 0.0398
    ##   ..$ corrected_effect   : num 0.115
    ##   ..$ corrected_effect_se: num 0.0535
    ##   ..$ test_difference    : num -2.36
    ##   ..$ p_difference       : num 0.0184
    ##  $ LDSC               :List of 11
    ##   ..$ h2_exp           : num 0.244
    ##   ..$ h2_exp_se        : num 0.0107
    ##   ..$ int_exp          : num 1.02
    ##   ..$ h2_out           : num 0.14
    ##   ..$ h2_out_se        : num 0.00975
    ##   ..$ int_out          : num 1.01
    ##   ..$ gcov             : num 0.0355
    ##   ..$ gcov_se          : num 0.00709
    ##   ..$ rg               : num 0.192
    ##   ..$ int_crosstrait   : num -0.0011
    ##   ..$ int_crosstrait_se: num 0.0063
    ##  $ GeneticArchitecture:List of 2
    ##   ..$ polygenicity       : num 0.00701
    ##   ..$ perSNP_heritability: num 3.02e-05

``` r
# MR + correction
unlist(A[["MRcorrection"]])
```

    ##     observed_effect  observed_effect_se    corrected_effect corrected_effect_se 
    ##          0.08560606          0.03980045          0.11457176          0.05353010 
    ##     test_difference        p_difference 
    ##         -2.35838564          0.01835461

  - **Example B**

<!-- end list -->

``` r
# observed effect
B[["MRcorrection"]]$observed_effect
```

    ## [1] 0.251392

``` r
# corrected effect
B[["MRcorrection"]]$corrected_effect
```

    ## [1] 0.2358383

``` r
# difference p-value
B[["MRcorrection"]]$p_difference
```

    ## [1] 1.070587e-05

``` r
# LDSC results
unlist(B[["LDSC"]])
```

    ##            h2_exp         h2_exp_se           int_exp            h2_out 
    ##       0.397445344       0.049379931       0.998226893       0.006929744 
    ##         h2_out_se           int_out              gcov           gcov_se 
    ##       0.023324039       0.994038341       0.073534339       0.020728556 
    ##                rg    int_crosstrait int_crosstrait_se 
    ##       1.401176698       0.361450785       0.006300000

## Runtime

Example A \~ 4 minutes 30 seconds

Example B \~ 3 minutes 40 seconds

The runtime can be influenced by the size of the summary statistics
files, the approach used for pruning (distance vs LD) but also by `s`,
the number of simulations used for the sampling strategy to estimate the
variance of the corrected causal effect and the covariance between
observed and corrected effects (type `?MRlap` to get details). However,
reducing the value of `s` is strongly discouraged as it can lead to an
innacurate estimation of the corrected effect standard error.

<font color="grey"><small> Results from analyses performed on a MacBook
Pro (Early 2015) - Processor : 2.9 GHz Intel Core i5 - Memory : 8 GB
1867 MHz DDR3.</font> </small>

## Citation

<!--- If you use the `MRlap` package, please cite:

[Ninon Mounier, Zoltán Kutalik, bGWAS: an R package to perform Bayesian Genome Wide Association Studies, Bioinformatics](https://doi.org/10.1093/bioinformatics/btaa549) --->

## Contact

<mounier.ninon@gmail.com>
