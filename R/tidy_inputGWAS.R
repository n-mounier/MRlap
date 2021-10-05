###### Function to create a nice data.frame for any type of input GWAS ######



# #' Tidy input GWAS
# #'
# #' From the GWAS arguments (exposure/outcome) of the main MRlap function,
# #' create a nice/tidy data.frame that can be used by all other functions.
# #'
# #' @param GWAS xx
# #'
# #' @inheritParams MRlap
# #' @export
# NOT EXPORTED



tidy_inputGWAS <- function(GWAS, K_t=NA, P_t=NA, verbose=FALSE){

  if(verbose) cat("# Preparation of the data... \n")

  GWASnames = list(SNPID = c("rsid", "snpid", "snp", "rnpid", "rs"),
                   CHR = c("chr"),
                   POS = c("pos"),
                   ALT = c("a1", "alt", "alts"),
                   REF = c("ref", "a0", "a2"),
                   BETA = c("beta", "b", "beta1", "or"),
                   SE = c("se", "std"),
                   Z = c("z", "zscore"),
                   N = c("n", "neff"),
                   Ncases = c("n_cases", "ncases", "n_case", "ncase"),
                   Ncontrols = c("n_controls", "ncontrols", "n_control", "ncontrol"))


  if(is.character(GWAS)) {

    # First, does the file exists ?
    if(!file.exists(GWAS)) stop("the file does not exist", call. = FALSE)
    if(verbose) cat(paste0("The file used as input is: \"",
                           strsplit(GWAS, "/")[[1]][length(strsplit(GWAS, "/")[[1]])], "\".  \n"))

    # Check colnames...
    HeaderGWAS = colnames(data.table::fread(GWAS, nrows = 1, showProgress = FALSE, data.table=F))


  } else if(is.data.frame(GWAS)){ # if data.frame
      # we want a data.frame (tidyr), not a data.table
    if(data.table::is.data.table(GWAS)) GWAS = as.data.frame(GWAS)
    if(verbose) cat(paste0("The data.frame used as input is: \"",
                           attributes(GWAS)$GName, "\".  \n"))
    HeaderGWAS = colnames(GWAS)

  }


  if(verbose & !is.na(K_t)) cat("   case-control data, analysis will be perform on the liability scale\n")


  HeaderGWAS = tolower(HeaderGWAS)


  if(all(!HeaderGWAS %in% GWASnames[["SNPID"]])) stop("no SNPID column", call. = FALSE)
  if(sum(HeaderGWAS %in% GWASnames[["SNPID"]])>1) stop("multiple SNPID columns, please provide only one", call. = FALSE)

  tmp = paste0("   SNPID column, ok")

  if(all(!HeaderGWAS %in% GWASnames[["CHR"]])) stop("no CHR column", call. = FALSE)
  tmp = c(tmp, "CHR column, ok")

  if(all(!HeaderGWAS %in% GWASnames[["POS"]])) stop("no POS column", call. = FALSE)
  tmp = c(tmp, "POS column, ok")

  if(all(!HeaderGWAS %in% GWASnames[["ALT"]])) stop("no ALT column", call. = FALSE)
  if(sum(HeaderGWAS %in% GWASnames[["ALT"]])>1) stop("multiple ALT columns, please provide only one", call. = FALSE)
  tmp = c(tmp, "ALT column, ok")

  if(all(!HeaderGWAS %in% GWASnames[["REF"]])) stop("no REF column", call. = FALSE)
  if(sum(HeaderGWAS %in% GWASnames[["REF"]])>1) stop("multiple REF columns, please provide only one", call. = FALSE)
  tmp = c(tmp, "REF column, ok")


  if(all(!HeaderGWAS %in% GWASnames[["Z"]])){
    if(!all(!HeaderGWAS %in% GWASnames[["BETA"]]) & !all(!HeaderGWAS %in% GWASnames[["SE"]])){
      tmp = c(tmp, "BETA column, ok")
      tmp = c(tmp, "SE column, ok")
      getZ=TRUE
    } else {
      stop("no effect (BETA/SE or Z) column(s)", call. = FALSE)
    }
  } else if(!all(!HeaderGWAS %in% GWASnames[["Z"]])){
    tmp = c(tmp, "Z column, ok")
    getZ=FALSE
  } else {
    stop("no effect (BETA/SE or Z) column(s)", call. = FALSE)
  }
  if(sum(HeaderGWAS %in% GWASnames[["BETA"]])>1) stop("multiple BETA columns, please provide only one", call. = FALSE)
  if(sum(HeaderGWAS %in% GWASnames[["SE"]])>1) stop("multiple SE columns, please provide only one", call. = FALSE)
  if(sum(HeaderGWAS %in% GWASnames[["Z"]])>1) stop("multiple Z columns, please provide only one", call. = FALSE)


  # sample size : if linear (or case-control without prevalences), should be N or Neff
  if(is.na(K_t)){
    if(all(!HeaderGWAS %in% GWASnames[["N"]])) stop("no N column", call. = FALSE)
    tmp = c(tmp, paste0("N column, ok \n"))
    if(sum(HeaderGWAS %in% GWASnames[["N"]])>1) stop("multiple N columns, please provide only one", call. = FALSE)
  } else {
  #               if case-control with prevalence, should be N or Neff and/or (Ncases + Ncontrols)
    getNeff = TRUE
    getNtot = TRUE
    if(all(HeaderGWAS %in% GWASnames[["N"]])){
      tmp = c(tmp, paste0("N column, ok \n"))
      getNeff=FALSE
    }
    if(!all(!HeaderGWAS %in% GWASnames[["Ncases"]]) & !all(!HeaderGWAS %in% GWASnames[["Ncontrols"]])){
      tmp = c(tmp, paste0("N (cases and controls) columns, ok \n"))
      getNtot=FALSE
    }
  if(getNeff & getNtot) stop("no N column", call. = FALSE)
  }
    # if(sum(HeaderGWAS %in% GWASnames[["N"]]) + sum(HeaderGWAS %in% GWASnames[["Ncases"]])>1) stop("multiple N columns (effective sample size + case / control sample sizes), please provide only one", call. = FALSE)
    # if(sum(HeaderGWAS %in% GWASnames[["N"]]) + sum(HeaderGWAS %in% GWASnames[["Ncontrols"]])>1) stop("multiple N columns (effective sample size + case / control sample sizes), please provide only one", call. = FALSE)


  if(verbose) cat(paste(tmp, collapse= " - "))





  # if headers ok, and file or data.frame get GWASData
  if(is.character(GWAS)) {
    # Get the full data
    GWASData = tibble::as_tibble(data.table::fread(GWAS, showProgress = FALSE, data.table=F))
    attributes(GWASData)$GName = basename(GWAS)

  } else if(is.data.frame(GWAS)){ # if data.frame
    # add attribute GName to the data.frame, to be re-used in other subfunctions
    GWASData=tibble::as_tibble(GWAS)
    rm(GWAS)
  }

  # use col numbers because of different lower/upper possibilities
  SNPID = match(HeaderGWAS, GWASnames[["SNPID"]])
  SNPID = which(!is.na(SNPID))[1]
  CHR = match(HeaderGWAS, GWASnames[["CHR"]])
  CHR = which(!is.na(CHR))[1]
  POS = match(HeaderGWAS, GWASnames[["POS"]])
  POS = which(!is.na(POS))[1]
  ALT = match(HeaderGWAS, GWASnames[["ALT"]])
  ALT = which(!is.na(ALT))[1]
  REF = match(HeaderGWAS, GWASnames[["REF"]])
  REF = which(!is.na(REF))[1]
  BETA = match(HeaderGWAS, GWASnames[["BETA"]])
  BETA = which(!is.na(BETA))[1]
  SE = match(HeaderGWAS, GWASnames[["SE"]])
  SE = which(!is.na(SE))[1]
  ZSTAT = match(HeaderGWAS, GWASnames[["Z"]])
  ZSTAT = which(!is.na(ZSTAT))[1]
  # for the sample size, if case-control, make sure to have Neff and Ntot
  if(is.na(K_t)){
    N = match(HeaderGWAS, GWASnames[["N"]])
    N = which(!is.na(N))[1]

    colNumbers = c(SNPID, CHR, POS, ALT, REF, BETA, SE, ZSTAT, N)
    colNames = c("rsid", "chr", "pos", "alt", "ref", "beta", "se", "Z", "N")
    colNames = colNames[!is.na(colNumbers)]
    colNumbers = colNumbers[!is.na(colNumbers)]


    GWASData %>%
      dplyr::select(dplyr::all_of(colNumbers)) %>%
      stats::setNames(colNames) -> GWASData_clean

  } else {
    if(!getNeff & getNtot){
      N = match(HeaderGWAS, GWASnames[["N"]])
      N = which(!is.na(N))[1]

      colNumbers = c(SNPID, CHR, POS, ALT, REF, BETA, SE, ZSTAT, N)
      colNames = c("rsid", "chr", "pos", "alt", "ref", "beta", "se", "Z", "N")
      colNames = colNames[!is.na(colNumbers)]
      colNumbers = colNumbers[!is.na(colNumbers)]


      GWASData %>%
        dplyr::select(dplyr::all_of(colNumbers)) %>%
        stats::setNames(colNames) %>%
        dplyr::mutate(Ntot = .data$N/(4 * P_t * (1-P_t))) -> GWASData_clean
    } else if(getNeff & !getNtot){
      N_cases = match(HeaderGWAS, GWASnames[["Ncases"]])
      N_cases = which(!is.na(N_cases))[1]
      N_controls = match(HeaderGWAS, GWASnames[["Ncontrols"]])
      N_controls = which(!is.na(N_controls))[1]

      colNumbers = c(SNPID, CHR, POS, ALT, REF, BETA, SE, ZSTAT, N_cases, N_controls)
      colNames = c("rsid", "chr", "pos", "alt", "ref", "beta", "se", "Z", "Ncases", "Ncontrols")
      colNames = colNames[!is.na(colNumbers)]
      colNumbers = colNumbers[!is.na(colNumbers)]


      GWASData %>%
        dplyr::select(dplyr::all_of(colNumbers)) %>%
        stats::setNames(colNames) %>%
        dplyr::mutate(Ntot = .data$Ncases+.data$Ncontrols,
                      N = 4*.data$Ncases*.data$Ncontrols/.data$Ntot) -> GWASData_clean
    } else if(!getNeff & !getNtot){
      N = match(HeaderGWAS, GWASnames[["N"]])
      N = which(!is.na(N))[1]
      N_cases = match(HeaderGWAS, GWASnames[["Ncases"]])
      N_cases = which(!is.na(N_cases))[1]
      N_controls = match(HeaderGWAS, GWASnames[["Ncontrols"]])
      N_controls = which(!is.na(N_controls))[1]

      colNumbers = c(SNPID, CHR, POS, ALT, REF, BETA, SE, ZSTAT, N, N_cases, N_controls)
      colNames = c("rsid", "chr", "pos", "alt", "ref", "beta", "se", "Z", "N", "Ncases", "Ncontrols")
      colNames = colNames[!is.na(colNumbers)]
      colNumbers = colNumbers[!is.na(colNumbers)]


      GWASData %>%
        dplyr::select(dplyr::all_of(colNumbers)) %>%
        stats::setNames(colNames) %>%
        dplyr::mutate(Ntot = .data$Ncases+.data$Ncontrols) -> GWASData_clean
    }
  }






  #if no Z, but BETA and SE, calculate Z
  if(getZ){
    # if OR, use log(OR)/SE
    if("or" %in% HeaderGWAS){
      dplyr::mutate(Z = log(.data$beta)/.data$se,
                    beta=NULL, se=NULL) -> GWASData_clean
    } else  {
      GWASData_clean %>%
        dplyr::mutate(Z = .data$beta/.data$se,
                      beta=NULL, se=NULL) -> GWASData_clean
    }
  } else {
    GWASData_clean %>%
      dplyr::mutate(beta=NULL, se=NULL) -> GWASData_clean
  }

  # Get standardised effects for MR
  # + LDSC needs a p-value column
  if(!is.na(K_t)){

    delta_t = ((K_t^2 * (1-K_t)^2) / (P_t * (1-P_t))) *
      ( 1 / ( dnorm( qnorm(1-K_t)) )^2 )

    i_t = dnorm( qnorm(1-K_t)) / K_t
    t= qnorm(1-K_t)
    theta_t = ( i_t * (P_t-K_t) / (1-K_t )) * (i_t * (P_t-K_t) / (1-K_t ) - t )


    GWASData_clean %>%
      dplyr::mutate(std_beta = sqrt(delta_t)*.data$Z/sqrt(.data$Ntot+delta_t*theta_t*.data$Z^2),
                    std_SE = 1/.data$Z * std_beta,
                    #Z = std_beta/std_SE,
                    p = 2*stats::pnorm(-abs(.data$Z))) -> GWASData_clean
  } else {
    GWASData_clean %>%
      dplyr::mutate(std_beta = .data$Z/sqrt(.data$N),
                    std_SE = 1/sqrt(.data$N),
                    p = 2*stats::pnorm(-abs(.data$Z))) -> GWASData_clean
  }
  res=GWASData_clean
  return(res)

}
