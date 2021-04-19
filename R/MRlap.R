#' MRlap - main function
#'
#' Performs cross-trait LD score regression, IVW-MR analysis and provide a correction
#' that simultaneously accounts for biases due to the overlap between the exposure and
#' outcome samples, the use of weak instruments and Winner’s curse.
#'
#'
#' @param exposure The path to the file containing the GWAS summary statistics for the exposure,
#'        or a \code{data.frame} (character, or \code{data.frame})
#' @param exposure_name The name of the exposure trait, \code{default="exposure"} (character)
#' @param outcome  The path to the file containing the GWAS summary statistics for the exposure,
#'        or a \code{data.frame} (character, or \code{data.frame})
#' @param outcome_name The name of the outcome trait, \code{default="outcome"} (character)
#' @param ld The path to the folder in which the LD scores used in the analysis are located.
#'        Expects LD scores formated as required by the original LD score regression software.  (character)
#' @param hm3 The path to a file of SNPs with alt, ref alleles and rsid used to allign alleles across traits
#'        (character)
#' @param MR_threshold The threshold used to select strong instruments for MR, should be lower
#'        than 1e-5, \code{default=5e-8} (numeric)
#' @param MR_pruning_dist The distance used for pruning MR instruments (in Kb), should be between 10 and 1000,
#'        \code{default=500} (numeric)
#' @param MR_pruning_LD The LD threshold (r2) used for pruning MR instruments, should be between 0 and 1
#'        (if 0, distance-based pruning is used), \code{default=0} (numeric)
#' @param MR_reverse The p-value used to exclude MR instruments that are more strongly associated with the outcome
#'        than with the exposure,\code{default=1e-3} (numeric)
#' @param s The number of simulations used for sampling strategy to estimate the variance of the corrected causal
#'        effect and the covariance between observed and corrected effects \code{default=10,000} (numeric)
#' @param save_logfiles  A logical indicating if log files from LDSC should be saved,
#'        \code{default=FALSE}
#' @param verbose  A logical indicating if information on progress should be reported,
#'        \code{default=TRUE}

#' @details
#' \code{exposure} and \code{outcome} are required arguments.
#' The input file / data.frame should contain the following
#' columns (lower or upper case) : \cr
#' SNPID (rs numbers) should be : \code{rs}, \code{rsid}, \code{snp}, \code{snpid}, \code{rnpid} \cr
#' CHR (chromosome) should be :  \code{chr} \cr
#' POS (position) should be :  \code{pos} \cr
#' ALT (effect allele) should be : \code{a1}, \code{alt}, \code{alts} \cr
#' REF (reference allele) should be : \code{a2}, \code{a0}, \code{ref} \cr
#' Z (z-score) should be : \code{Z}, \code{zscore} \cr
#' N (sample size) should be :  \code{N} \cr
#' If Z is not present, it can be calculated from BETA and SE. \cr
#' BETA should be : \code{b}, \code{beta}, \code{beta1} \cr
#' SE should be : \code{se}, \code{std} \cr
#'
#' Setting \code{s} to a smaller value is strongly discouraged, it can lead to an innacurate
#' estimatation of the corrected effect SE, therefore affecting the results of the difference testing
#' between observed and corrected effects.
#'
#' @importFrom rlang .data
#' @export

MRlap <- function(exposure,
                  exposure_name = NULL,
                  outcome,
                  outcome_name = NULL,
                  ld,
                  hm3,
                  MR_threshold = 5e-8,
                  MR_pruning_dist = 500,
                  MR_pruning_LD = 0,
                  MR_reverse = 1e-3,
                  s=10000,
                  save_logfiles = FALSE,
                  verbose = TRUE) {


  # Path where the analysis has been launched
  InitPath = getwd()
  on.exit(setwd(InitPath))

  StartTime =  proc.time()

  # platform identification (not needed)
  platform = .Platform$OS.type

  ## Be chatty ?
  if(!is.logical(verbose)) stop("verbose : should be logical", call. = FALSE)
  if(!is.logical(save_logfiles)) stop("save_logfiles : should be logical", call. = FALSE)


  if(verbose) cat("<<< Preparation of analysis >>> \n")

  ### check the parameters ###
  if(verbose) cat("> Checking parameters \n")


  if(is.data.frame(exposure)){# add attribute GName to the data.frame, to be re-used in other subfunctions
    attributes(exposure)$GName =  deparse(substitute(exposure)) # get the "name" of the object used as an argument in the function
  }
  if(is.data.frame(outcome)){# add attribute GName to the data.frame, to be re-used in other subfunctions
    attributes(outcome)$GName =  deparse(substitute(outcome)) # get the "name" of the object used as an argument in the function
  }


  ## ld_wd + hm3
  if (is.character(ld)){
    if(!dir.exists(ld)) stop("ld : the folder does not exist", call. = FALSE)
    # get absolute path
    ld = normalizePath(ld)
  } else stop("ld : wrong format, should be character", call. = FALSE)
  if (is.character(hm3)){
    if(!file.exists(hm3)) stop("hm3 : the file does not exist", call. = FALSE)
    # get absolute path
    hm3 = normalizePath(hm3)
  } else stop("hm3 : wrong format, should be character", call. = FALSE)


  ## MR_threshold -> should not be larger than 10-5, can only be more stringent
  if(!is.numeric(MR_threshold)) stop("MR_threshold : non-numeric argument", call. = FALSE)
  if(MR_threshold>10^-5) stop("MR_threshold : superior to the threshold limit", call. = FALSE)

  if(verbose) cat("The p-value threshold used for selecting MR instruments is:", format(MR_threshold, scientific = T), "\n")

  ## MR_pruning_dist
  if(!is.numeric(MR_pruning_dist)) stop("MR_pruning_dist : non-numeric argument", call. = FALSE)
  if(MR_pruning_dist<10) stop("MR_pruning_dist : should be higher than 10Kb", call. = FALSE)
  if(MR_pruning_dist>1000) stop("MR_pruning_dist : should be lower than 1Mb", call. = FALSE)


  if(verbose) cat("The distance used for pruning MR instruments is: ", MR_pruning_dist, "Kb \n")

  ## MR_pruning_LD
  if(!is.numeric(MR_pruning_LD)) stop("MR_pruning_LD : non-numeric argument", call. = FALSE)
  if(MR_pruning_LD<0) stop("MR_pruning_LD : should be positive", call. = FALSE)
  if(MR_pruning_LD>1) stop("MR_pruning_LD : should not be larger than 1", call. = FALSE)

  if(MR_pruning_LD>0){
    if(verbose) cat("The LD threshold used for pruning MR instruments is:", MR_pruning_LD, "\n")
  } else {
    if(verbose) cat("Distance-based pruning will be used for MR instruments \n")
  }

  # 0 : Tidy input GWAS
  if(!is.null(exposure_name) & !is.character(exposure_name)) stop("exposure_name : non-character argument", call. = FALSE)
  if(!is.null(outcome_name) & !is.character(outcome_name)) stop("outcome_name : non-character argument", call. = FALSE)

  if(verbose) cat(paste0("> Processing exposure ", ifelse(is.null(exposure_name), "", paste0("(",exposure_name,") ")) ,"summary statistics... \n"))
  if(is.null(exposure_name)) exposure_name="exposure"
  exposure_data = tidy_inputGWAS(exposure, verbose)
  if(verbose) cat(paste0("> Processing outcome ", ifelse(is.null(outcome_name), "", paste0("(",outcome_name,") ")) ,"summary statistics... \n"))
  if(is.null(outcome_name)) outcome_name="outcome"
  outcome_data = tidy_inputGWAS(outcome, verbose)



  ### 1 : run cross-trait LDSC ###
  if(verbose) cat("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> \n")
  if(verbose) cat("<<< Performing cross-trait LDSC >>>  \n")
  # returns h2 exposure - SE h2 exposure - cross-trait intercept - SE cross-trait intercept
  LDSC_results = run_LDSC(exposure_data, exposure_name, outcome_data, outcome_name, ld, hm3, save_logfiles, verbose)
  # -> h2_LDSC, h2_LDSC_se, lambda, lambda_se (for correction)
  #     int_exp, int_out,  h2_out, h2_out_se, rgcov, rgcov_se, rg

  # 2 : run IVW-MR
  if(verbose) cat("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> \n")
  if(verbose) cat("<<< Running IVW-MR >>>  \n")
  # returns alpha - SE alpha - instruments (needed for corrected effect SE)
  MR_results = run_MR(exposure_data, outcome_data, MR_threshold,
                      MR_pruning_dist, MR_pruning_LD, MR_reverse,
                      verbose)
  # -> alpha_obs, alpha_obs_se, n_exp, n_out, IVs
  # 3 : get corrected effect
  if(verbose) cat("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> \n")
  if(verbose) cat("<<< Estimating corrected effect >>>  \n")
  correction_results = with(c(MR_results, LDSC_results),
    get_correction(IVs, lambda, lambda_se, h2_LDSC, h2_LDSC_se,
                                      alpha_obs, alpha_obs_se,
                                      n_exp, n_out, MR_threshold, verbose, s))
  # -> alpha_corrected, alpha_corrected_se, cov_obs_corrected, test_diff, p_diff
  #    pi_x, sigma2_x


  tmp = "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n"


  ### we're done! ###
  Time = as.integer((proc.time()-StartTime)[3])
  minutes <- as.integer(trunc(Time/60))
  seconds <- Time - minutes * 60
  if(verbose) cat("Runtime of the analysis: ", minutes, " minute(s) and ", seconds, " second(s).  \n")


  # results -> list of 3
  # [[1]] "MR correction"
  # [[2]] "LDsc" : h2X / seh2X / h2Y / seh2Y /
  # [[3]] "GeneticArchitecture" : pi / sigma
  results_MR=with(c(MR_results, correction_results),
               list("observed_effect" = alpha_obs,
                    "observed_effect_se" = alpha_obs_se,
                    "corrected_effect" = alpha_corrected,
                    "corrected_effect_se" = alpha_corrected_se,
                    "test_difference" = test_diff,
                    "p_difference" = p_diff))

  results_LDSC=with(LDSC_results,
                  list("h2_exp" = h2_LDSC ,
                       "h2_exp_se" = h2_LDSC_se,
                       "int_exp" = int_exp,
                       "h2_out" = h2_out,
                       "h2_out_se" = h2_out_se,
                       "int_out" = int_out,
                       "gcov" = rgcov,
                       "gcov_se" = rgcov_se,
                       "rg" = rg,
                       "int_crosstrait" = lambda,
                       "int_crosstrait_se" = lambda_se))

  results_GeneticArchitecture=with(correction_results,
                    list("polygenicity" = pi_x,
                         "perSNP_heritability" = sigma2_x))

  results = list(MRcorrection = results_MR,
                 LDSC = results_LDSC,
                 GeneticArchitecture = results_GeneticArchitecture)
  return(results)
}
