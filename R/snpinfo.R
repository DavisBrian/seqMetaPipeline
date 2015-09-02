# This function provides a basic check that the SNPInfo file is a suitable for
# seqMeta.  It check that the required column names are in the SNPInfo data
# frame and that they are of the correct type. 
#
checkSNPInfo <- function(snpinfo, snpNames=NULL, aggregateBy=NULL, filterBy=NULL, chrName=NULL) {
  if (is.null(snpNames) & is.null(aggregateBy)) {
    #check that there are no factor 
    if (any(sapply(snpinfo, is.factor))) {
      stop("Factor columns not allowed in snpinfo.  Convert to characters.")
    }
  }
  
  # check snpNames
  if (!is.null(snpNames)) {
    if (!(snpNames %in% colnames(snpinfo))) {
      stop(paste0("snpNames: ", snpNames, " is not a column name of the provided SNPInfo file."))
    }
    if (typeof(snpinfo[ , snpNames]) != "character") {
      stop(paste0("snpNames: ", snpNames, "is of type ", typeof(snpinfo[ , snpNames]), ".  Must be of type character"))
    }
  }
  
  # check aggregateBy
  if (!is.null(aggregateBy)) {
    if (!(aggregateBy %in% colnames(snpinfo))) {
      stop(paste0("aggregateBy: ", aggregateBy, " is not a column name of the provided SNPInfo file."))
    }
    if (typeof(snpinfo[ , aggregateBy]) != "character") {
      stop(paste0("aggregateBy: ", aggregateBy, "is of type ", typeof(snpinfo[ , aggregateBy]), ".  Must be of type character"))
    }
  }
  
  # check filterBy
  if (!is.null(filterBy)) {
    if (any(!(filterBy %in% colnames(snpinfo)))) {
      stop(paste0("filterBy: ", filterBy[!(filterBy %in% colnames(snpinfo))], " is not a column name of the provided SNPInfo file."))
    }
    if (length(filterBy) == 1L) {
      if (typeof(snpinfo[ , filterBy]) != "logical") {
        stop(paste0("filterBy: ", filterBy, "is of type ", typeof(snpinfo[ , chrName]), ".  Must be of type logical"))
      }      
    } else if (length(filterBy) >= 1L) {
      test <- sapply(snpinfo[, snpinfo.filterBy], typeof)
      if (any(test != "logical")) {
        stop(paste0("filterBy: ", names(test)[test != "logical"], "is of type ", test[test != "logical"], ".  Must be of type logical"))
      }      
    }
  }  
  
  #check chrName
  if (!is.null(chrName)) {
    if (!(chrName %in% colnames(snpinfo))) {
      stop(paste0("chrName: ", chrName, " is not a column name of the provided SNPInfo file."))
    }
    if (typeof(snpinfo[ , chrName]) != "character") {
      stop(paste0("chrName: ", chrName, "is of type ", typeof(snpinfo[ , chrName]), ".  Must be of type character"))
    }
  }  
  
  return(invisible(NULL))   
}

