#' MRlap - main function
#'
#' Performs a cross-trait LD score regression, IVW-MR analysis and
#' CORRECTION!
#'
#'
#' @param exposure The path to the file containing the GWAS summary statistics for the exposure,
#'        or a \code{data.frame} (character, or \code{data.frame})
#' @param outcome  The path to the file containing the GWAS summary statistics for the exposure,
#'        or a \code{data.frame} (character, or \code{data.frame})
#' @param ld The path to the folder in which the LD scores used in the analysis are located.
#'        Expects LD scores formated as required by the original LD score regression software.  (character)
#' @param wld The path to the folder in which the LD scores weights used in the analysis are located.
#'        Expects LD scores formated as required by the original LD score regression software.  (character)
#' @param MR_threshold The threshold used to select strong instruments for MR, should be lower
#'        than 1e-5, \code{default=5e-8} (numeric)
#' @param MR_pruning_dist The distance used for pruning MR instruments (in Kb), should be between 10 and 1000,
#'        \code{default=500} (numeric)
#' @param MR_pruning_LD The LD threshold used for pruning MR instruments, should be between 0 and 1
#'        (if 0, distance-based pruning is used), \code{default=0} (numeric)
#' @param verbose  A logical indicating if information on progress should be reported,
#'        \code{default=TRUE}

#' @details
#' \code{exposure} and \code{outcome} are required arguments.
#' The input file / data.frame should contain the following
#' columns : \cr
#' SNPID (rs numbers) should be : \code{rs}, \code{rsid}, \code{snp}, \code{snpid}, \code{rnpid} \cr
#' ALT (effect allele) should be : \code{a1}, \code{alt}, \code{alts} \cr
#' REF (reference allele) should be : \code{a2}, \code{a0}, \code{ref} \cr
#' Z (z-score) should be : \code{z}, \code{Z}, \code{zscore} \cr
#' N (sample size) should be :
#' If Z is not present, it can be calculated from BETA and SE. \cr
#' BETA should be : \code{b}, \code{beta}, \code{beta1} \cr
#' SE should be : \code{se}, \code{std} \cr

MRlap <- function(exposure,
                  outcome,
                  ld_wd,
                  MR_threshold = 1e-6,
                  MR_pruning_dist = 500,
                  MR_pruning_LD = 0,
                  verbose = TRUE) {


  # Path where the analysis has been launched
  InitPath = getwd()
  on.exit(setwd(InitPath))

  StartTime =  proc.time()

  # platform identification
  # if platform is windows, an path should be provided for Z-matrices
  platform = .Platform$OS.type
  if(platform=="windows" && Z_matrices=="~/ZMatrices/") stop("Windows operating system, please provide a file path to the \"Z_matrices\" argument", call. = FALSE)


  ## Be chatty ?
  if(!is.logical(verbose)) stop("verbose : should be logical", call. = FALSE)


  # initialization of log_info file
  log_info = c()

  if(verbose) cat("<<< Preparation of analysis >>> \n")

  ### check the parameters ###
  if(verbose) cat("> Checking parameters \n")



  ## ld_wd
  if (is.character(Z_matrices)){

    if(!dir.exists(Z_matrices)) stop("Z_matrices : the folder do not exist", call. = FALSE)

    # get absolute path
    Z_matrices = normalizePath(Z_matrices)

    if(!file.exists(file.path(Z_matrices, "ZMatrix_Full.csv.gz"))) stop("No \"ZMatrix_Full.csv.gz\" file in specified Z_matrices folder", call. = FALSE)
    if(!file.exists(file.path(Z_matrices, "ZMatrix_MR.csv.gz"))) stop("No \"ZMatrix_MR.csv.gz\" file in specified Z_matrices folder", call. = FALSE)
    if(!file.exists(file.path(Z_matrices, "AvailableStudies.tsv"))) stop("No \"AvailableStudies.tsv\" file in specified Z_matrices folder", call. = FALSE)


  } else stop("Z_matrices : wrong format, should be character", call. = FALSE)


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
    if(verbose) cat("The LD threshold used for pruning MR instruments is:", MR_pruning_LD)
  } else {
    if(verbose) cat("Distance-based pruning will be used for MR instruments \n")
  }

  # Tidy input GWAS
  if(verbose) cat("> Processing exposure summary statistics... \n")
  exposure_data = tidy_inputGWAS(exposure, verbose)
  if(verbose) cat("> Processing outcome summary statistics... \n")
  outcome_data = tidy_inputGWAS(exposure, verbose)



  ### 1 : run cross-trait LDSC ###
  if(verbose) cat("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> \n",
                  "<<< Performing cross-trait LDSC >>>  \n")
  run_LDSC(exposure_data, outcome_data)

  # 2 : run IVW-MR
  if(verbose) cat("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> \n",
                  "<<< Running IVW-MR >>>  \n")

  # 3 : get corrected effect
  if(verbose) cat("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> \n",
                  "<<< Estimating corrected effect >>>  \n")


  results=list()
  return(results)
}
