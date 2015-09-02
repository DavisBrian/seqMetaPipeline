# This is a wrapper function around prepScores, prepScoresX from the seqMeta 
# package.  seqMeta asumes that the people in the data phenotype
# data are the same as in the genotype matrix.  With missing data this is not
# generally the case.  Thus one must subset to the intersecting subjects. Also
# there is no gaurentee that all snps in the genotype matrix match the snpinfo
# file.  This function also trys to split the snpinfo file and genotype matrix
# if there are a combination of snps in the X chromosome and non-X chromosome.
#
createSeqMetaObject <- function(Z, FORMULA, pheno, FAMILY="gaussian", MALE=NULL,
                                SNPInfo=NULL, snpNames="Name", 
                                aggregateBy="gene", chrName="CHROM", kins=NULL, 
                                sparse = TRUE, verbose = FALSE) {
  
  if (is.null(FAMILY)) {
    fam <- gaussian()
  } else if (FAMILY == "gaussian") {
    fam <- gaussian()
  } else if (FAMILY == "binomial") {
    fam <- binomial()
  } else if (FAMILY == "survival") {
    fam <- FAMILY
  } else {
    stop("Only gaussian and binomial valid family types")
  }
  
  checkSNPInfo(SNPInfo, snpNames = snpNames, aggregateBy = aggregateBy, 
               filterBy = NULL, chrName = chrName)
  checkGenotypeMatrix(Z)
  checkPhenotype(pheno, FORMULA, idCol = NULL, genderCol=MALE)
  
  snps <- intersect(colnames(Z), SNPInfo[ , snpNames])
  subjects <- intersect(rownames(Z), rownames(pheno))
  
  #basic checks
  if (length(subjects) == 0L) {
    stop("No subjects in common between the genotype matrix and the phenotype data frame.")
  }
  if (length(snps) == 0L) {
    stop("No snps in common between the genotype matrix and the snpinfo file.")
  }
  
  if (length(subjects) != nrow(pheno)) {
    wmsg <- paste("Only", length(subjects), "of", nrow(pheno), "subjects found.", sep = " ") 
    message(wmsg)
    ZnotP <- setdiff(rownames(Z), rownames(pheno))
    if (length(ZnotP > 0L)) {
      ZnotPmsg <- paste(ZnotP, "does not have phenotype information.", sep = ' ')
      message(ZnotPmsg)
    }
    PnotZ <- setdiff(rownames(pheno), rownames(Z))
    if (length(PnotZ > 0L)) {
      PnotZmsg <- paste(PnotZ, "does not have genotype information.", sep=' ')
      message(PnotZmsg)
    }
  }
  if (length(snps) != ncol(Z)) {
    message(paste("Only", length(snps), "of", ncol(Z), "SNPs found.", sep=" "))
    ZnotSI <- setdiff(colnames(Z), SNPInfo[ , snpNames])
    if (length(ZnotSI > 0L)) {
      ZnotSImsg <- paste(ZnotSI, "does not have an SNP Info.", sep=' ')
      message(ZnotSImsg)
    }
  } 
  
  
  SI <- SNPInfo[SNPInfo[ , snpNames] %in% snps, ]
  chroms <- unique(SI[,  chrName])
  
  if (FAMILY == "survival") {
    # prepeCOx
    ARGS <- list(Z=Z[subjects, snps], 
                 formula=as.formula(FORMULA),
                 SNPInfo=SI,
                 snpNames=snpNames,
                 aggregateBy=aggregateBy,
                 data=pheno[subjects,],
                 verbose=verbose)    
    obj <- do.call(prepCox, args=ARGS)    
  } else {
    if (!is.null(MALE) & ("X" %in% chroms)){
      if (length(chroms) == 1L) {
        ARGS <- list(Z=Z[subjects, snps], 
                     formula=as.formula(FORMULA),
                     male=pheno[subjects, MALE],
                     family=fam,
                     SNPInfo=SI,
                     snpNames=snpNames,
                     aggregateBy=aggregateBy,
                     kins=kins,
                     sparse=sparse,
                     data=pheno[subjects,],
                     verbose=verbose)
        obj <- do.call(prepScoresX, args=ARGS)
      } else {
        # run X
        X <- (SI[ , chrName] == "X")
        snps.X <- intersect(colnames(Z), SI[X , snpNames])
        ARGS <- list(Z=Z[subjects, snps.X], 
                     formula=as.formula(FORMULA),
                     male=pheno[subjects, MALE],
                     family=fam,
                     SNPInfo=SI[X, ],
                     snpNames=snpNames,
                     aggregateBy=aggregateBy,
                     kins=kins,
                     sparse=sparse,
                     data=pheno[subjects,],
                     verbose=verbose)            
        obj.X <- do.call(prepScoresX, args=ARGS)
        
        # run NotX
        snps.NotX <- intersect(colnames(Z), SI[!X , snpNames])
        ARGS <- list(Z=Z[subjects, snps.NotX], 
                     formula=as.formula(FORMULA),
                     family=fam,
                     SNPInfo=SI[!X, ],
                     snpNames=snpNames,
                     aggregateBy=aggregateBy,
                     kins=kins,
                     sparse=sparse,
                     data=pheno[subjects,],
                     verbose=verbose)      
        obj.NotX <- do.call(prepScores, args=ARGS)      
        
        # put back together
        obj <- do.call(seqMeta:::c.seqMeta, args=list(obj.NotX, obj.X))
      }   
    } else {      
      ARGS <- list(Z=Z[subjects, snps], 
                   formula=as.formula(FORMULA),
                   family=fam,
                   SNPInfo=SI,
                   snpNames=snpNames,
                   aggregateBy=aggregateBy,
                   kins=kins,
                   sparse=sparse,
                   data=pheno[subjects, ],
                   verbose=verbose) 
      obj <- do.call(prepScores, args=ARGS)   
    }     
  }
  
  return(obj)  
}

# This is a wrapper function around burdenMeta which provides a bit more error
# checking for common failure modes.  It also adds in the estimated MAC to the
# results.
# TBD: Add in the calculated MAC from the actual data
#
runBurdenMeta <- function(cohort, snpinfo=NULL, wts=1, snpNames="Name", 
                          aggregateBy="gene", mafRange=c(0,0.5), 
                          verbose=FALSE, chrName=NULL) {
  
  #  check we have deqMeta or skatMeta cohort objects
  if (class(cohort) != "seqMeta" & class(cohort) != "skatCohort") {
    stop("cohort is not a seqMeta object!")
  }
  
  checkSNPInfo(snpinfo, snpNames=snpNames, aggregateBy=aggregateBy, 
               filterBy=NULL, chrName=chrName)
  
  MAC.calc <- CalculateMAC(cohort, snpinfo=snpinfo, snpNames=snpNames, 
                           aggregateBy=aggregateBy, mafRange=mafRange)
  #  MAC.calc[is.na(MAC.calc)] <- 0
  MAC <- data.frame(MAC.calc)
  MAC$gene <- rownames(MAC)
  
  if (!is.null(chrName)) {
    CHR <- unique(snpinfo[, c(chrName, aggregateBy)])
    MAC <- merge(MAC, CHR, by.x="gene", by.y=aggregateBy, all.x=TRUE, all.y=FALSE, sort=FALSE)    
  }
  
  # get N for each gene.
  N <- ldply(cohort, .fun=function(x) {x$n})
  colnames(N) <- c("gene", "N")
  
  res <- burdenMeta(cohort, SNPInfo=snpinfo, wts=wts, snpNames=snpNames, 
                    aggregateBy=aggregateBy, mafRange=mafRange, verbose=verbose)
  results <- merge(res, N, by="gene", all.x=TRUE, all.y=FALSE, sort=FALSE)
  results$MAC.est <- 2*results$N*results$cmafUsed
  results <- merge(results, MAC, by="gene", all.x=TRUE, all.y=FALSE, sort=FALSE)
  
  return(results)  
}

# This is a wrapper function around skatMeta which provides a bit more error
# checking for common failure modes.  It also adds in the estimated MAC to the
# results.
# TBD: Add in the calculated MAC from the actual data
#
runSkatMeta <- function(cohort, snpinfo=NULL, wts=function(maf){dbeta(maf,1,25)}, 
                        method="saddlepoint", snpNames="Name", 
                        aggregateBy="gene", mafRange=c(0,0.5), verbose=FALSE, 
                        chrName=NULL) {
  
  #  check we have seqMeta or skatMeta cohort objects
  if (class(cohort) != "seqMeta" & class(cohort) != "skatCohort") {
    stop("cohort is not a seqMeta object!")
  }
  
  checkSNPInfo(snpinfo, snpNames=snpNames, aggregateBy=aggregateBy, 
               filterBy=NULL, chrName=chrName)
  
  MAC.calc <- CalculateMAC(cohort, snpinfo=snpinfo, snpNames=snpNames, 
                           aggregateBy=aggregateBy, mafRange=mafRange)
  MAC <- data.frame(MAC.calc)
  MAC$gene <- rownames(MAC)
  if (!is.null(chrName)) {
    CHR <- unique(snpinfo[, c(chrName, aggregateBy)])
    MAC <- merge(MAC, CHR, by.x="gene", by.y=aggregateBy, all.x=TRUE, all.y=FALSE, sort=FALSE)    
  }
  
  # get N for each gene.
  N <- ldply(cohort, .fun=function(x) {x$n})
  colnames(N) <- c("gene", "N")
  
  # this is how SKAT should be run.  However, there is a bug in seqMeta 1.3
  # where if all the snps are filtered out by mafRange it will generate an
  # error and crash.  The solution is to include all snps in the SKAT test
  # but assign a weight of 0 to all snps outside the mafRange.  The downside
  # is that the number of SNPs and the CMAF outputted are not correct using
  # the work-around but you can obtain the correct number of SNPs/CMAF by 
  # combining with the burden test with the same snpinfo and MAF threshold.
  #     SKAT <- skatMeta(cohort, SNPInfo=snpinfo, snpNames=snpNames, 
  #                      aggregateBy=aggregateBy, mafRange=c(0,0.05), verbose=verbose)
  res <- skatMeta(cohort, SNPInfo=snpinfo, wts=wts, method=method, 
                  snpNames=snpNames, aggregateBy=aggregateBy, 
                  verbose=verbose) 
  results <- merge(res, N, by="gene", all.x=TRUE, all.y=FALSE, sort=FALSE)
  results$MAC.est <- 2*results$N*results$cmaf
  results <- merge(results, MAC, by="gene", all.x=TRUE, all.y=FALSE, sort=FALSE)    
  
  return(results)    
}



# This is a wrapper function around singlesnpMeta which provides a bit more error
# checking for common failure modes.  It also adds in the estimated MAC to the
# results.
# TBD: Add in the calculated MAC from the actual data
#
runSingleSNP <- function(cohort, snpinfo=NULL, snpNames="Name", 
                         aggregateBy="gene", studyBetas=TRUE,verbose=FALSE) {
  
  #  check we have deqMeta or skatMeta cohort objects
  if (class(cohort) != "seqMeta" & class(cohort) != "skatCohort") {
    stop("cohort is not a seqMeta object!")
  }
  
  checkSNPInfo(snpinfo, snpNames=snpNames, aggregateBy=aggregateBy, 
               filterBy=NULL, chrName=NULL)
  
  # get N for each gene.
  N <- ldply(cohort, .fun=function(x) {x$n})
  colnames(N) <- c("gene", "N")
  
  res <- singlesnpMeta(cohort, SNPInfo=snpinfo, snpNames=snpNames, 
                       aggregateBy=aggregateBy, studyBetas=studyBetas, 
                       verbose=verbose)
  results <- merge(res, N, by="gene", all.x=TRUE, all.y=FALSE, sort=FALSE)
  results$MAC.est <- 2*results$N*results$maf
  if ("MAC" %in% colnames(snpinfo)) {
    results <- merge(results, snpinfo[, c(snpNames, "MAC")], by.x="Name", 
                     by.y=snpNames, all.x=TRUE, all.y=FALSE, sort=FALSE) 
  } else {
    results$MAC <- NA
  }  
  return(results)     
}
