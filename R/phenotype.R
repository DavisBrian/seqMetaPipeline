# This function provides a basic check that the phenotype file is a suitable
# for seqMeta.  It checks that the required column names are in the phenotype
# data frame. 
#
checkPhenotype <- function(p, pformula, idCol=NULL, genderCol=NULL) {
  if (!is.null(idCol)) {
    if (!(idCol %in% colnames(p))) {
      stop("ID column not in the phenotype data frame.")
    }
    if (any(duplicated(p[ , idCol]))) {
      stop("Duplicated phenotype ids.")
    }
  }
  
  if(!is.null(genderCol)) {
    if (!(genderCol %in% colnames(p))) {
      msg <- paste(genderCol, "not found in phenotype file.", sep=" ")
      stop(msg)
    }
    gtype <- typeof(p[ , genderCol])
    g <- unique(p[ , genderCol])
    if(gtype == "integer") {
      if (!all(g %in% c(0, 1))) {
        stop("Gender must be (0/1 or F/T) indicating female/male.")
      }      
    } else if(gtype == "character") {
      if (!all(g %in% c("F", "M"))) {
        stop("Gender must be (0/1 or F/T) indicating female/male.")
      }        
    }     
  } else {
    wmsg <- "No column given to identify Males in phenotype file.  No special 
    handling of the X chromosme."
    warning(wmsg)
  }
  
  tf <- terms(as.formula(pformula), data=p)
  if (grepl("^Surv\\(", pformula)) {
    cn <- attr(tf, "term.labels")
  } else {
    cn <- rownames(attr(tf, "factors"))
  }
  
  if (any(!(cn %in% colnames(p)))) {
    msg <- paste("Formula varaibles:", cn[!(cn %in% colnames(p))], "not found in phenotype file.", sep=" ")
    stop(msg)
  } 
  return(invisible(NULL)) 
}

# This function reduces a data set to only the variables used in a model
# removing subjects with missing data.  Also, it makes the row names of
# the resulting data fram the subject identifier
#
# p: a data frame containing the variables in the model
#
# formula: a character vector which can be coered to an object of class 
#          "formula" (with as.formula): a symbolic description of the model to
#          be fitted. The details of model specification are given under 
#          'Details' in the "lm" help file.
#
# id: (optional) colunm name identifier of the subjects
#
# gender: (optional) colunm name identifier for the gender classification of
#         the subjects.
#
# returns: data frame with only the columns specified in the formula and with
#          the (optional) row names as the subject identifier.
#     
reducePheno <- function(p, pformula, id=NULL, gender=NULL) {
  checkPhenotype(p, pformula, idCol=id, genderCol=gender)   
  
  if (!is.null(id)) {
    if (!(id %in% colnames(p))) {
      stop("ID column not in the phenotype data frame.")
    } else {
      rownames(p) <- p[ , id]
      p <- p[, -match(id, colnames(p))]
    }
  } 
  
  if (!is.null(gender)) {
    gtype <- typeof(p[ , gender])
    g <- unique(p[ , gender])
    if(gtype == "character") {
      #     if (all.equal(g, c("F", "M"))) {
      if (all.equal(sort(g), c("F", "M"))) {
        MF <- p[, gender] == "M"
        p[, gender] <- MF
      } else {
        stop("Unable to convert gender.  Gender must be (0/1 or F/T) indicating female/male.")
      }       
    }                                   
  }
  
  #   tf <- terms(as.formula(pformula), data=p)
  #   cn <- rownames(attr(tf, "factors"))
  tf <- terms(as.formula(pformula), data=p)
  if (grepl("^Surv\\(", pformula)) {
    tcn <- attr(tf, "term.labels")
    fcn <- min(match(tcn, colnames(p)), na.rm=TRUE)
    cn <- union(colnames(p)[1:(fcn-1)], tcn)    
  } else {
    cn <- rownames(attr(tf, "factors"))
  }
  
  
  if ((!is.null(gender)) & (!(gender %in% cn))) {
    cn <- c(cn, gender)    
  }
  #   if (!(gender %in% cn)) {
  #     cn <- c(cn, gender)
  #   }
  
  if (any(!(cn %in% colnames(p)))) {
    msg <- paste("Formula varaibles:", cn[!(cn %in% colnames(p))], "not found in phenotype file.", sep=" ")
    stop(msg)
  } else {
    if (!is.null(cn)) {
      pheno <- na.omit(p[, cn])
    } else {
      pheno <- na.omit(p)
    }    
  }
  
  return(pheno)
}

summarize_phenotype <- function(pheno, gender=NULL) {
  if(!is.null(gender)) {
    gtype <- typeof(pheno[ , gender])
    g <- unique(pheno[ , gender])
    if(gtype == "integer") {
      if (!all(g %in% c(0, 1))) {
        stop("Gender must be (0/1 or F/T) indicating female/male.")
      }  
      
      MALES <- which(pheno[ , gender] == 1L)
    } else if(gtype == "character") {
      if (!all(g %in% c("F", "M"))) {
        stop("Gender must be (0/1 or F/T) indicating female/male.")
      } 
      MALES <- which(pheno[ , gender] == "M")    
    } else if(gtype == "logical") {
      MALES <- which(pheno[ , gender] == TRUE)
    } else {
      wmsg("Gender must be (0/1 or F/T) indicating female/male.  Cannot calculate gener specific phenotype summaries")
      warning(wmsg)
      MALES <- NULL
    }
  } else {
    wmsg("No column given to identify Males in phenotype file.  Cannot calculate gener specific phenotype summaries")
    warning(wmsg)
    MALES <- NULL
  }
  
  pooled.summary <- CalculatePhenotypeSummary(pheno)
  if (is.null(MALES)) {
    male.summary <- NA
    female.summary <- NA
  } else {
    male.summary <- CalculatePhenotypeSummary(pheno[MALES,])
    female.summary <- CalculatePhenotypeSummary(pheno[-MALES,])    
  }
  
  return(list(Pooled=pooled.summary, Male=male.summary, Female=female.summary)) 
}

CalculatePhenotypeSummary <- function(p) {
  
  cols <- sapply(p, function(x) { (is.numeric(x) | is.integer(x))})
  
  N <- colSums(!is.na(p[, cols]))
  Min <- sapply(p[, cols], min, na.rm=TRUE)
  Q1 <- sapply(p[, cols], quantile, probs=0.25, na.rm=TRUE)
  Median <- sapply(p[, cols], median, na.rm=TRUE)
  Mean <- colMeans(p[, cols], na.rm=TRUE)
  SD <- sapply(p[, cols], sd, na.rm=TRUE)
  Q3 <- sapply(p[, cols], quantile, probs=0.75, na.rm=TRUE)
  Max <- sapply(p[, cols], max, na.rm=TRUE)
  
  p.summary <- rbind(N, Min, Q1, Median, Mean, SD, Q3, Max)
  
  return(p.summary)
}
