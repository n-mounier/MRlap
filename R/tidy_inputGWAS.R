###### Function to create a nice data.frame for any type of input GWAS######



# #' Tidy input GWAS
# #'
# #' From the GWAS arguments (exposure/outcime) of the main MRlap function,
# #' create a nice/tidy data.frame that can be used by all other functions.
# #'
# #' @inheritParams MRlap
# NOT EXPORTED



tidy_inputGWAS <- function(GWAS, verbose=FALSE){

  if(verbose) cat("# Preparation of the data... \n")


  GWASnames = list(SNPID = c("rsid", "snpid", "snp", "rnpid", "rs"),
                   ALT = c("a1", "alt", "alts"),
                   REF = c("ref", "a0", "a2"),
                   BETA = c("beta", "b", "beta1"),
                   SE = c("se", "std"),
                   Z = c("z", "zscore"),
                   N = c("n", "N"))


  if(is.character(GWAS)) {

    # First, does the file exists ?
    if(!file.exists(GWAS)) stop("the file does not exist", call. = FALSE)
    if(verbose) cat(paste0("The file used as input is: \"",
                           strsplit(GWAS, "/")[[1]][length(strsplit(GWAS, "/")[[1]])], "\".  \n"))

    # Check colnames...
    HeaderGWAS = colnames(data.table::fread(GWAS, nrows = 1, showProgress = FALSE, data.table=F))


  } else if(is.data.frame(GWAS)){ # if data.frame
    # add attribute GName to the data.frame, to be re-used in other subfunctions
    attributes(GWAS)$GName =  deparse(substitute(GWAS)) # get the "name" of the object used as an argument in the function
    # we want a data.frame (tidyr), not a data.table
    if(data.table::is.data.table(GWAS)) GWAS = as.data.frame(GWAS)
    if(verbose) cat(paste0("The data.frame used as input is: \"",
                           attributes(GWAS)$GName, "\".  \n"))
    HeaderGWAS = colnames(GWAS)

  }



  HeaderGWAS = tolower(HeaderGWAS)


  if(all(!HeaderGWAS %in% GWASnames[["SNPID"]])) stop("no SNPID column", call. = FALSE)
  # how to deal with multiple rsid / snpid columns ???
  # here, we don't care, we need at least one
  tmp = paste0("   SNPID column, ok")

  if(all(!HeaderGWAS %in% GWASnames[["ALT"]])) stop("no ALT column", call. = FALSE)
  tmp = c(tmp, paste0("ALT column, ok"))

  if(all(!HeaderGWAS %in% GWASnames[["REF"]])) stop("no REF column", call. = FALSE)
  tmp = c(tmp, paste0("REF column, ok"))


  # if beta + se
  if(!all(!HeaderGWAS %in% GWASnames[["BETA"]]) & !all(!HeaderGWAS %in% GWASnames[["SE"]])){
    if(!all(!HeaderGWAS %in% GWASnames[["Z"]])){
      tmp = c(tmp, paste0("Z column, ok"))
      getZ=FALSE
    } else { # by default, use Z, only use beta/se if z not provided
      getZ=TRUE
    }

    tmp = c(tmp, paste0("BETA column, ok"))
    tmp = c(tmp, paste0("SE column, ok"))
  } else if(!all(!HeaderGWAS %in% GWASnames[["Z"]])){
    tmp = c(tmp, paste0("Z column, ok"))
  } else {
    stop("no effect (BETA/SE or Z) column(s)", call. = FALSE)
  }

  if(all(!HeaderGWAS %in% GWASnames[["N"]])) stop("no N column", call. = FALSE)
  tmp = c(tmp, paste0("N column, ok"))


  if(verbose) cat(paste(tmp, collapse= " - "))





  # if headers ok, and file or data.frame get GWASData
  if(is.character(GWAS)) {
    # Get the full data
    GWASData = as_tibble(data.table::fread(GWAS, showProgress = FALSE, data.table=F))
    attributes(GWASData)$GName = basename(GWAS)

  } else if(is.data.frame(GWAS)){ # if data.frame
    # add attribute GName to the data.frame, to be re-used in other subfunctions
    GWASData=as_tibble(GWAS)
    rm(GWAS)
  }

  # by default, always use the first column of input GWAS matching the names
  # use col numbers because of different lower/upper possibilities
  SNPID = match(HeaderGWAS, GWASnames[["SNPID"]])
  SNPID = which(!is.na(SNPID))[1]
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
  N = match(HeaderGWAS, GWASnames[["N"]])
  N = which(!is.na(N))[1]


  colNumbers = c(SNPID, ALT, REF, BETA, SE, ZSTAT, N)
  colNames = c("rsid", "alt", "ref", "beta", "se", "z", "N")
  colNames = colNames[!is.na(colNumbers)]
  colNumbers = colNumbers[!is.na(colNumbers)]


  GWASData %>%
    select(all_of(colNumbers)) %>%
    stats::setNames(colNames) -> GWASData_clean

  #if no Z, but BETA and SE, calculate Z
  if(getZ){
    GWASData_clean %>%
      mutate(z_obs = .data$beta/.data$se) -> GWASData_clean
  }

  GWASData_clean %>%
    mutate(std_beta = Z/sqrt(N)) -> GWASData_clean

  res=GWASData_clean
  return(res)

}
