checkGenotypeMatrix <- function(Z) {
  if (!is.matrix(Z)) {
    stop("Genotype must be in a numeric matrix.")
  }
  
  # make sure the data is numeric or integer.
  if ((typeof(Z) != "integer") & (typeof(Z) != "double")) {
    stop("Genotype matrix must be numeric (integer or double).")    
  }
  
  return(invisible(NULL)) 
  
}
