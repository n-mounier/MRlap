# MRlap 0.0.3.1 (2023-10-11)

## Changes

- instrument selection:   

User can use local installation of Plink for LD clumping with a custom reference panel (using `MR_plink` and `MR_bfile` options)    
User can turn of automatic pruning and instead provide own list of variants to use as Instrumental Variables (using `do_pruning` and `user_SNPsToKeep` options)    

- input data format   

the chromosome and position columns are not needed anymore if LD pruning is performed using a local version of PLINK (still required for distance pruning, and for LD pruning using the API, to do it chromosome by chromosome and work with smaller subsets of SNPs).    
As a consequence, the HLA region is not excluded by default before performing MR (if you prefer to exclude it, please do it before launching MRlap)


# MRlap 0.0.3.0 (2022-11-18)

## Bug fixes
- Parametric bootstrap when the exposure has relatively low heritability (and/or when the corresponding standare error is relatively large)   

In such a case, negative heritability values could be generated in the pametric bootstrap. We are now ensuring that this doesn't lead to any error in the bootstrap procedure and added a warning message when this happens, as very low heritability could  compromise the results.

## Changes
- rsids  

The rsids of the of instruments used are now reported.


# MRlap 0.0.2.0 (2021-10-26)

## Bug fixes
- Excluding IVs more strongly associated with the outcome than with the exposure   

There was a typo leading to an underestimation of the t-statistics, and IVs more strongly associated with the outcome than with the exposure were not correctly identified.

## Changes
- rsids  

The rsids of the of instruments used are now reported.

# MRlap 0.0.1.0 (2021-10-15)

## Changes
- Case-control GWAS    

It is now possible to use case-control GWAS for the exposure and/or the outcome datasets.
If (at least) one of the datasets is coming from a case-control GWAS, the Sample size column should correspond to the total sample size (not the effective sample size).   
The analysis is performed on the observed scale (explaining why heritability estimates will differ from the ones usually obtained on the liability scale) and the approach is the same as the one used for continuous traits.    

- optimal number of simulations to estimate SE/COV

The number of simulation used to estimate the standard error of the corrected effect and the covariance between corrected and observed effect is no longer a parameter that needs to be chosen before running the analysis. This could be a problem if the value was too small, the SE/COV estimates were too noisy to correcly test for the difference between observed and corrected effects. The optimal value is now automatically determined when running the analysis (using the coefficients of variation for the variance of the corrected effect and covariance between corrected and observed effect). 

- p-value estimation for observed and corrected effects / number of IVs

The p-values are now reported for both observed and corrected effects, as well as the number of instruments used.

# MRlap 0.0.0.9 (2021-03-04)

initial version of the package, as described in the preprint

<!--- 
## Bug fixes

## New functions

## Documentation

## Error messages

## Performance


--->  
