is_try_error <- function(x) inherits(x, "try-error")

# This Calculates the MAC without imputing missing values from the cohort
# object and the phenotype data used.  This is done by first calling
# AddMACtoSNPInfo and adding the MAC from the current model to the snpinfo file
#
CalculateMAC <- function(cohort, snpinfo, snpNames, aggregateBy, mafRange) {
  #  check we have deqMeta or skatMeta cohort objects
  if (class(cohort) != "seqMeta" & class(cohort) != "skatCohort") {
    stop("cohort is not a seqMeta object!")
  }
  # check we have the correct columns in our snpinfo file
  if (!(snpNames %in% colnames(snpinfo))) {
    msg <- paste("snapNames: ", snpNames, " does not match any column name in SNPInfo", sep='')
    stop(msg)
  }
  if (!(aggregateBy %in% colnames(snpinfo))) {
    msg <- paste("aggregateBy: ", aggregateBy, " does not match any column name in SNPInfo", sep='')
    stop(msg)
  }
  
  genes <- names(cohort)
  MAC <- integer(length(genes))
  names(MAC) <- genes
  MAC[]<-NA
  
  if ("MAC" %in% colnames(snpinfo)){
    for (gene in genes) {
      x<-cohort[[gene]]
      si <- snpinfo[(snpinfo[, aggregateBy] %in% gene), c(snpNames, aggregateBy, "MAC")]
      
      snps <- unique(names(x$maf)[(x$maf > min(mafRange)) & (x$maf <= max(mafRange))])
      if (length(snps) > 0L) {
        idx <- match(snps, si[, snpNames])
        MAC[gene] <- sum(si[idx, "MAC"], na.rm=TRUE)
      }  
    }
  }  
  
  return(MAC)
  
}

# This will calculate the MAC based on the current genotype matrix and 
# phenotype file.  This is later used by the burden tests to calcualte MAC
# for each gene
#
AddMACtoSNPInfo <- function(Z, pheno, snpinfo, snpNames = "Name") { 
  if (!(snpNames %in% colnames(snpinfo))) {
    msg <- paste("snpNames: ", snpNames, " does not match any column name in SNPInfo", sep='')
    stop(msg)
  }
  
  subjects <- intersect(rownames(Z), rownames(pheno))
  if (length(subjects) == 0L) {
    warning("No subjects in common between the genotype matrix and the phenotype data frame.")
  }
  
  MAC <- data.frame(colnames(Z[subjects, ]), colSums(Z[subjects, ], na.rm=TRUE), stringsAsFactors=FALSE)
  colnames(MAC) <- c("SNP", "MAC")
  
  si <- merge(snpinfo, MAC, by.x=snpNames, by.y="SNP", all.x=TRUE, all.y=FALSE, sort=FALSE)
  return(si)
}

