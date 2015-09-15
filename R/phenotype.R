# phenotype

#' @title Create a phenotype object
#'   
#' @description This function creates a phenotype object.
#' 
#' @param data a data frame (or object coercible by as.data.frame to a data 
#'   frame) containing the variables in the model.
#' @param formula an object of class "formula" (or one that can be coerced to 
#'   that class): a symbolic description of the model to be fitted.The details 
#'   of model specification are given under 'Details' in the "lm" help file.
#' @param id (optional) column name identifier of the unique subjects in
#'   \code{data}.  If given the phenotype sample ids will become the rownames
#'   of \code{data}. If \code{id} is \code{NULL} the suject id's are assumed
#'   to be the rownames of \code{data}.
#' @param gender (optional) column name identifier for the gender in \code{data}
#' @param include (optional) character vector of the subjects in \code{id} to be
#'   included in the analysis.  See Details.
#' @param exclude (optional) character vector of the subjects in \code{id} to be
#'   excluded in the analysis.  See Details.
#' @param reduce logical.  Should the dataset be reduced to only columns used in
#' formula.  See Details
#'   
#' @details \code{data} and \code{formula} are similar to what is needed for
#'   \code{lm}.  If a formula is not specifed, no further analysis will run.
#'   
#'   If the \code{id} is not specified it is assumed that the rownames in
#'   \code{data} are the unique subject identifier.  If \code{id} is specifed
#'   then \code{rownames(data)} will be set to \code{id}.  Thus the subject id's
#'   must be unique and without duplication.
#'   
#'   If the gender column is not specified it will be set to \code{NULL}. 
#'   
#'   Typically either \code{include} or \code{exclude}, or both, is set to
#'   \code{NULL}. It is important to note that return \code{include} and
#'   \code{exclude} are taken from the subjects in \code{data} and not from the input parameters.
#'   [TBD: explain include/exclude]
#'   
#'   If \code{reduce == TRUE} then \code{data} is reduced to a data.frame
#'   containing the variables used in \code{formula} plus \code{gender}. See
#'   get_all_vars for specifics.
#'   
#'  @seealso get_all_vars 
#'   
#' @return an object of class 'phenotype'.  This is a list typicaly used to feed
#'   into an analysis function.
#'   
#'   An object of class 'phenotype' is a list containing at least the following
#'   components:
#'   
#' \itemize{
#'  \item data data frame of the data for analysis 
#'  \item formula formula to be used in further analysis  
#'  \item gender column name containing the gender of the subjects
#'  \item include the subjects in \code{data} that will be included is further analysis.
#'  \item exclude the subjects in \code{data} that will be excluded is further analysis.    
#' }
#' 
#' It is important to note that return \code{include} and \code{exclude} are taken from data and not from the input parameters. 
#'   
#' @export
#
# [TBD]
#  -verbose option?
#  -add genderChar??? something to demote which character is "MALE/FEMALE"
#  -add "family" (gaussian/binomial/survival)
#  -add print method to show the meta data
#  -add "problems" (a la readr)
phenotype <- function(data, formula, family, id, gender, 
                      include=NULL, exclude=NULL) {
  
  # check data
  if (!is.data.frame(data)) {
    stop("Error in parameter 'data': Must be a data.frame or a class which extends a data.frame.")
  }
  
  # check formula
  if (missing(formula)) {
    stop("'formula' must be specified!")
  }
  # - make sure all varaibles are in the data frame
  vars <- all.vars(as.formula(formula))
  vars <- vars[!(vars %in% ".")]
  if (all(vars %in% colnames(data))) {
    mf <- get_all_vars(formula, data)
    if (nrow(na.omit(mf)) <= 1L) {
      stop("'formula' results in empty 'data'")
    }
  } else {
    stop("One or more formula variables not in data")
  }

  
  # check family
  if (missing(family) || is.null(family) || length(family) == 0L || is.na(family)) {
    stop("'family' must be specified!")
  }
  if (is.character(family)) {
    if (length(family) > 1L) {
      stop("Family can be only one of: 'gaussian', 'binomial', or 'cox'")
    }
    if (!all(family %in% c('gaussian', 'binomial', 'cox'))) {
      stop("Only family 'gaussian', 'binomial', or 'cox' currently supported")
    }
  } else {
    stop("'family' must be of type character.")
  }

  # check id
  if (missing(id)) {
    warning("id not specified.  Using row names")
    id <- ".id"
    data$.id <- rownames(data)
  } else {
    if (length(id) != 1L) {
      stop("Only one id column can be specified.")
    } else {
      if (!(id %in% colnames(data))) {
        stop(paste0(id, " is not a column in data."))
      }
    }
    if (!is.character(data[ , id])) {
      warning("Converting id column to character.")
      data[ , id] <- as.character(data[ , id]) 
    }
    if (anyDuplicated(data[ , id])) {
      stop("Duplicate phenotype ids found.")
    }
    if (anyNA(data[ , id])) {
      stop("Missing phenotype ids")
    }
  }

  # check gender
  if (missing(gender)) {
    gender <- NULL
  } else {
    if (length(gender) != 1L) {
      stop("Only one gender column can be specified.")
    } else {
      if (!(gender %in% colnames(data))) {
        stop(paste0(gender, " is not a column in data."))
      }
    }
    if (anyNA(data[ , gender])) {
      stop("Missing gender data")
    }
    # [TBD]  - is a type than can be grouped
    # [TBD] - has no more than 2 groups
    # [TBD] - can be converted to TRUE/FALSE??    
  }

  
  subjects_all <- data[ , id]
  
  # check include
  if (!is.null(include)) {
    if (!is.character(include)) {
      stop("include must be a character vector")
    }
    if (!all(include %in% data[ ,id])) {
      warning("Not all ids in 'include' are in 'data'")
    }
    if (!all(include %in% data[ ,id])) {
      warning("Not all ids in 'include' are in 'data'")
    }
    subjects_include <- intersect(include, data[ , id])    
    data <- data[(data[, id] %in% subjects), , drop = FALSE]
  }

  # check exclude
  if (!is.null(exclude)) {
    if (!is.character(exclude)) {
      stop("exclude must be a character vector")
    }
    if (!all(exclude %in% data[ ,id])) {
      warning("Not all ids in 'exclude' are in 'data'")
    }
    subjects_exclude <- setdiff(data[ , id], exclude)
    data <- data[(data[, id] %in% subjects), , drop = FALSE]
  }
  
  # check null model can be built with the parameters given.
  
#   nullmodel <- if (family == "binomial" || family == "gaussian") {
#     try(glm(formula = as.formula(formula), family = family, data = data[ , -idCol]))
#   } else if (family == "cox") {
#     try(coxph(formula = as.formula(formula), data = data[ , -idCol]))
#   } else {
#     stop("Unknown family type.  Only 'gaussian', 'binomial', and 'cox' are currently supported.")
#   }
#   if (is_try_error(nullmodel)) {
#     stop("nullmodel can not be built.")
#   }
  
#   data <- na.omit(data[ , cols])
#   
#   subjects_include <- data[ , id]
#   
#   subjects_exclude <- setdiff(subjects_all, subjects_include)
#   if (length(subjects_exclude) == 0L) {
#     subjects_exclude <- NULL
#   }
#   
  new_class <- class(data)
  
#   structure(
#     data,
#     formula = formula,
#     family = family,
#     idCol = id,
#     genderCol = gender,
#     included = subjects_include,
#     excluded = subjects_exclude,
#     class = unique(c("phenotype", new_class))
#   )
    structure(
      data,
      formula = formula,
      family = family,
      idCol = id,
      genderCol = gender,
      class = unique(c("phenotype", new_class))
    )
  
}

#' @rdname phenotype
#' @export
is_phenotype <- function(x) inherits(x, "phenotype")

# get functions  ---------------------------------------------------------------

#' @export
get_subjects.phenotype <- function(x, excluded = FALSE) {
  # [TBD] check exlude is  a logical of length 1
  stopifnot(length(excluded) == 1L)
  if (excluded) {
    attr(x, "excluded")
  } else {
    intersect(x[[attr(x, "idCol")]], attr(x, "included"))
  }
}

#' @export
get_formula <- function(x) { attr(x, "formula") }

#' @export
get_family <- function(x) { attr(x, "family") }

#' @export
get_idCol <- function(x) { attr(x, "idCol") }

#' @export
get_genderCol <- function(x) { attr(x, "genderCol") }

#' @export
get_included <- function(x) { attr(x, "included") }

#' @export
get_excluded <- function(x) { attr(x, "excluded") }


# single verbs -----------------------------------------------------------------

#' @export
reduce.phenotype <- function(p, common) {
  
  dropped <- setdiff(get_subjects(p), common$subjects)
  if (length(dropped) == 0L) {
    attr(p, "dropped") <- NA
  } else {
    attr(p, "dropped") <- dropped
  }
  
  o <- match(common$subjects, p[ , get_idCol(p)])
  p[o, ]
}

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
