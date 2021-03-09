#' Exposure
#'
#' Subset of the original dataset containing the estimated effect of SNPs on the exposure
#'
#' @format A data frame with 750,000 rows and 13 variables:
#' \describe{
#'   \item{chr}{chromosome}
#'   \item{rsid}{rsid of the SNP}
#'   \item{pos}{position}
#'   \item{ref}{reference allele for the SNP}
#'   \item{alt}{effect allele for the SNP}
#'   \item{af}{allele frequency}
#'   \item{info}{imputation quality}
#'   \item{beta}{estimated effect size for the SNP}
#'   \item{se}{standard error of the estimated effect size for the SNP}
#'   \item{z}{z-score for the SNP}
#'   \item{minuslog10p}{-log10(p) for the SNP}
#'   \item{p}{p-value for the SNP}
#'   \item{N}{sample size}
#' }
#' @source \url{}
"SmallExposure_Data"


#' Outcome
#'
#' Subset of the original dataset containing the estimated effect of SNPs on the outcome
#'
#' @format A data frame with 750,000 rows and 13 variables:
#' \describe{
#'   \item{chr}{chromosome}
#'   \item{rsid}{rsid of the SNP}
#'   \item{pos}{position}
#'   \item{ref}{reference allele for the SNP}
#'   \item{alt}{effect allele for the SNP}
#'   \item{af}{allele frequency}
#'   \item{info}{imputation quality}
#'   \item{beta}{estimated effect size for the SNP}
#'   \item{se}{standard error of the estimated effect size for the SNP}
#'   \item{z}{z-score for the SNP}
#'   \item{minuslog10p}{-log10(p) for the SNP}
#'   \item{p}{p-value for the SNP}
#'   \item{N}{sample size}
#' }
#' @source \url{}
"SmallOutcome_Data"





