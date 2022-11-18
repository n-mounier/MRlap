
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
:information\_source: `MRlap` has been updated to version 0.0.3. Note
that in previous versions the exclusion of IVs more strong associated
with the outcome than with the exposure was not correcltly performed.  
Check the [NEWS](NEWS.md) to learn more about what has been modified!

## Overview

`MRlap` is an R-package to perform two-sample Mendelian Randomisation
(MR) analyses using (potentially) overlapping samples, relying only on
GWAS summary statistics. MR estimates can be subject to different types
of biases due to the overlap between the exposure and outcome samples,
the use of weak instruments and winner’s curse. Our approach
simultaneously accounts and corrects for all these biases, using
cross-trait LD-score regression (LDSC) to approximate the overlap.
Estimating the corrected effect using our approach can be performed as a
sensitivity analysis: if the corrected effect do not significantly
differ from the observed effect, then IVW-MR estimate can be safely
used. However, when there is a significant difference, corrected effects
should be preferred as they should be less biased, independently of the
sample overlap.  
Note that we are working with standardised effects. This means that the
causal effect estimates are in units of Standard Deviation (SD). The
causal effect estimates correspond to the SD change in the outcome for
one SD increase in the exposure.  
It is possible to use case-control GWASs, but it is important to note
that our method assumes that sample overlap should be independend of
case-control status. Morever, the analysis needs to be performed on the
observed scale. This means that GWAS results (from a linear model using
case-control status or from a logisitic regression) should be provided
alongside the total (number of cases + number of controls) sample size.
In such cases, the heritability estimates reported in the results will
be different to the ones that are usually estimated (on the liability
scale), this is normal. The causal effect estimates correspond to the SD
changes / increases on the observed scale.  
This package builds up on the
[`GenomicSEM`](https://github.com/GenomicSEM/GenomicSEM/) R-package to
perform cross-trait LDSC and the
[`TwoSampleMR`](https://github.com/MRCIEU/TwoSampleMR/) R-package for
inverse-variance weighted (IVW-)MR analysis (and instruments pruning).

There is only one function available:

-   **`MRlap()`**  
    main function that performs LDSC, IVW-MR analysis and provides a
    corrected causal effect estimate.

More details about its usage can be found in the
[manual](doc/MRlap-manual.pdf).

## Installation

You can install the current version of `MRlap` with:

``` r
# Directly install the package from github
# install.packages("remotes")
remotes::install_github("n-mounier/MRlap", )
library(MRlap)
```

<!--- Note: using remotes instead of devtools leads to re-build the package
and apparently, it may be a problem with R 3.4 and macOS, 
see https://stackoverflow.com/questions/43595457/alternate-compiler-for-installing-r-packages-clang-error-unsupported-option/43943631#43943631 --->

## Usage

To run the analysis with `MRlap` different inputs are needed:

#### 1. The exposure and outcome GWAS summary statistics (`exposure` & `outcome`):

Can be a regular (space/tab/comma-separated) file or a gzipped file
(.gz) or a `data.frame`. Must contain the following columns, which can
have alternative names (not case sensitive):  
SNP-identifier: `rs` or `rsid`, `snp`, `snpid`, `rnpid`  
Alternate (effect) allele: `a1` or `alt`, `alts`  
Reference allele: `a2`, `a0`, `ref`  
Z-statistics: `Z`, `zscore`  
Sample size: `N`, `Neff`

If the Z-statistics is not present, it will be automatically calculated
from effect size and standard error, in which case the following columns
should be provided:  
Effect-size: `b`, `beta`, `beta1` , `or`  
Standard error: `se`, `std`

*If (at least) one of the datasets is coming from a case-control
GWAS:*  
… the Sample size column should correspond to the total sample size (not
the effective sample size!!).

#### 2. The input files for LDSC (`ld` & `hm3`):

These are needed by the
[`GenomicSEM`](https://github.com/GenomicSEM/GenomicSEM/) R-package.

-   ld:

> Expects LD scores formated as required by the original LD score
> regression software. Weights for the european population can be
> obtained by downloading the eur\_w\_ld\_chr folder in the link below
> (Note that these are the same weights provided by the original
> developers of LDSC):
> <https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v>

-   hm3:

> This file can be obtained from
> <https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v>.

### Analysis

Before running the examples, please make sure to have downloaded the
input files for LDSC. You may also need to modify the `ld` & `hm3`
parameters to indicate the correct paths.

-   **Example A**

``` r
# Using ~100K samples for BMI/SBP, with 0% of sample overlap
# (only weak instrument bias and winner's curse)
# Note that here the overlap is known (since we generated the data) but the MRlap
# function works even the overlap is unkown (overlap is *not* a parameter of the function) 
# as it uses cross-trait LDSC to approximate it
# (1,150,000 SNPs - stored in gzipped files)
BMI <- system.file("data/", "BMI_Data.tsv.gz", package="MRlap")
SBP <- system.file("data/", "SBP_Data.tsv.gz", package="MRlap")

# MR instruments will be selected using default parameter (5e-8) and distance-pruned (500Kb),
# No file will be saved.
A = MRlap(exposure = BMI,
          exposure_name = "BMI_100Ksample",
          outcome = SBP,
          outcome_name = "SBP_100Ksample",
          ld = "~/eur_w_ld_chr",
          hm3 = "~/w_hm3.noMHC.snplist")
```

<details>
<summary>
Show log
</summary>

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
    ##      corrected effect: 0.115 ( 0.0536 ) 
    ##      covariance between observed and corrected effect: 0.00215   
    ##            10000 simulations were used to estimate the variance and the covariance.
    ##  > Testing difference between observed and corrected effect...  
    ##  Runtime of the analysis:  3  minute(s) and  11  second(s).
    ```

</details>

-   **Example B**

``` r
# Using simulated exposure/outcome data 
# standard settings scenario, with 100% of sample overlap
# Note that here the overlap is known (since we generated the data) but the MRlap
# function works even the overlap is unkown (overlap is *not* a parameter of the function) 
# as it uses cross-trait LDSC to approximate it
# (~750,000 SNPs - stored as data.frames)

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
<summary>
Show log
</summary>

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
    ##      209 IVs with p < 5e-10  
    ##      0 IVs excluded - more strongly associated with the outcome than with the exposure, p < 1e-03  
    ##     Pruning : distance :  500 Kb  - LD threshold :  0.05  
    ##      38 IVs left after pruning  
    ##  > Performing MR  
    ##      IVW-MR observed effect: 0.217 ( 0.0235 ) 
    ##  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
    ##  <<< Estimating corrected effect >>>   
    ##  > Estimating genetic architecture parameters...  
    ##  > Estimating corrected effect...  
    ##      corrected effect: 0.199 ( 0.0265 ) 
    ##      covariance between observed and corrected effect: 0.000625  
    ##            7000 simulations were used to estimate the variance and the covariance.
    ##  > Testing difference between observed and corrected effect...  
    ##  Runtime of the analysis:  3  minute(s) and  09  second(s).
    ```

</details>

### Results

**`MRlap()`** returns a named list containing the following results:

-   MRcorrection

“observed\_effect” : IVW-MR observed causal effect estimate,  
“observed\_effect\_se” : IVW-MR observed causal effect standard error,  
“m\_IVs” : number of IVs used,  
“IVs” : rsid of IVs used,  
“observed\_effect\_p” : IVW-MR observed causal effect p-value,  
“corrected\_effect” : corrected causal effect estimate,  
“corrected\_effect\_se” : corrected causal effect standard error,  
“corrected\_effect\_p” : corrected causal effect p-value,  
“test\_difference” : test statistic used to test for differences between
observed and corrected effects,  
“p\_difference” : p-value corresponding to the test statistic used to
test for differences between observed and corrected effects.

-   LDSC

“h2\_exp” : exposure heritability estimate,  
“h2\_exp\_se” : exposure heritability standard error,  
“int\_exp” : exposure intercept estimate,  
“h2\_out” : outcome heritability estimate,  
“h2\_out\_se” : outcome heritability standard error,  
“int\_out” : outcome intercept estimate,  
“gcov” : genetic covariance estimate,  
“gcov\_se” : genetic covariance standard error,  
“rg” : genetic correlation estimate,  
“int\_crosstrait” : cross-trait intercept estimate,  
“int\_crosstrait\_se”: cross-trait intercept standard error.

-   GeneticArchitecture

“polygenicity” : exposure polygenicity estimate,  
“perSNP\_heritability” : exposure per-SNP heritability estimate.

##### Aditionnaly, if `save_logfiles=TRUE`, LDSC log files are created in the current working directory :

-   **<exposure_name>.log** - exposure cleaning/munging log file  

-   **<outcome_name>.log** - outcome cleaning/munging log file  

-   **<exposure_name>.sumstats.gz\_<outcome_name>.sumstats.gzldsc.log** -
    cross-trait LDSC log file

-   **Example A**

``` r
# structure of the results
str(A)
```

    ## List of 3
    ##  $ MRcorrection       :List of 10
    ##   ..$ observed_effect    : num 0.0856
    ##   ..$ observed_effect_se : num 0.0398
    ##   ..$ m_IVs              : int 39
    ##   ..$ IVs                : chr [1:39] "rs684382" "rs2051086" "rs543874" "rs6728726" ...
    ##   ..$ observed_effect_p  : num 0.0315
    ##   ..$ corrected_effect   : num 0.115
    ##   ..$ corrected_effect_se: num 0.0536
    ##   ..$ corrected_effect_p : num 0.0324
    ##   ..$ test_difference    : num -2.36
    ##   ..$ p_difference       : num 0.0185
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
    ##   ..$ polygenicity       : num 0.00702
    ##   ..$ perSNP_heritability: num 3.02e-05

``` r
# MR + correction
unlist(A[["MRcorrection"]])
```

    ##      observed_effect   observed_effect_se                m_IVs 
    ## "0.0856060582715485" "0.0398004451909771"                 "39" 
    ##                 IVs1                 IVs2                 IVs3 
    ##           "rs684382"          "rs2051086"           "rs543874" 
    ##                 IVs4                 IVs5                 IVs6 
    ##          "rs6728726"           "rs713586"          "rs6731302" 
    ##                 IVs7                 IVs8                 IVs9 
    ##          "rs2861683"          "rs1567959"          "rs6792892" 
    ##                IVs10                IVs11                IVs12 
    ##          "rs1492014"          "rs1320903"         "rs13130484" 
    ##                IVs13                IVs14                IVs15 
    ##          "rs2307111"          "rs2410763"         "rs16867703" 
    ##                IVs16                IVs17                IVs18 
    ##          "rs2235569"         "rs12528998"         "rs12202969" 
    ##                IVs19                IVs20                IVs21 
    ##          "rs1167796"          "rs2245368"          "rs2299381" 
    ##                IVs22                IVs23                IVs24 
    ##          "rs2439823"           "rs900144"             "rs6265" 
    ##                IVs25                IVs26                IVs27 
    ##         "rs10838777"          "rs7132908"         "rs10146997" 
    ##                IVs28                IVs29                IVs30 
    ##         "rs16951304"         "rs10083738"          "rs4889606" 
    ##                IVs31                IVs32                IVs33 
    ##          "rs1421085"          "rs4783718"          "rs4239060" 
    ##                IVs34                IVs35                IVs36 
    ##          "rs1916295"         "rs11653498"           "rs571312" 
    ##                IVs37                IVs38                IVs39 
    ##         "rs11672660"          "rs3810291"           "rs400140" 
    ##    observed_effect_p     corrected_effect  corrected_effect_se 
    ## "0.0314855203412581"  "0.114571580921726"  "0.053561613723008" 
    ##   corrected_effect_p      test_difference         p_difference 
    ## "0.0324306952903651"  "-2.35617980011097" "0.0184639782420551"

``` r
# in this case, we observed that the corrected effects points towards an underestimation
# of the observed effect estimate obtained using IVW (because when there is no sample 
# overlap winner's curse and weak instrument bias will bias the estimate towards the null)

# LDSC results
unlist(A[["LDSC"]])
```

    ##            h2_exp         h2_exp_se           int_exp            h2_out 
    ##       0.243865854       0.010663133       1.019459090       0.140424951 
    ##         h2_out_se           int_out              gcov           gcov_se 
    ##       0.009745983       1.011189974       0.035547929       0.007089347 
    ##                rg    int_crosstrait int_crosstrait_se 
    ##       0.192095266      -0.001101472       0.006300000

-   **Example B**

``` r
# observed effect
B[["MRcorrection"]]$observed_effect
```

    ## [1] 0.2167283

``` r
# corrected effect
B[["MRcorrection"]]$corrected_effect
```

    ## [1] 0.1985984

``` r
# difference p-value
B[["MRcorrection"]]$p_difference
```

    ## [1] 6.200377e-10

``` r
# in this case, we observed that the the observed effect estimate obtained using IVW 
# is overestimated because of the sample overlap. The true causal effect used for 
# simulating the data is 0.2 (bias for corrected effect is 3.5 folds lower).

# Exposure genetic architecture (estimated to get corrected effects)
unlist(B[["GeneticArchitecture"]])
```

    ##        polygenicity perSNP_heritability 
    ##        0.0008382092        0.0004034253

## Runtime

Example A \~ 3 minutes 40 seconds

Example B \~ 3 minutes 35 seconds

The runtime can be influenced by the size of the summary statistics
files, the approach used for pruning (distance vs LD) but also by the
number of simulations used for the sampling strategy to estimate the
variance of the corrected causal effect and the covariance between
observed and corrected effects (optimal number of simulations is
automatically determined in the analysis).

<font color="grey"><small> Results from analyses performed on a MacBook
Pro (2020) - Processor : 2 GHz Quad-Core Intel Core i5 - Memory : 16 GB
3733 MHz LPDDR4X.</font> </small>

## Citation

<!--- If you use the `MRlap` package, please cite:

[Ninon Mounier, Zoltán Kutalik, bGWAS: an R package to perform Bayesian Genome Wide Association Studies, Bioinformatics](https://doi.org/10.1093/bioinformatics/btaa549) --->

## Contact

<mounier.ninon@gmail.com>
