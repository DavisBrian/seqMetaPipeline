combineResults <- function(files) {
  
  combined.results <- NULL 
  
  for (FILE in files) {
    result.name <- load(FILE)
    if (result.name != "results") {
      results <- get(result.name)
      rm(list=result.name)  
    }
    combined.results <- rbind(combined.results, results)
  }
  
  return(combined.results)
}


writeResults <- function(results, outdir, test, baseName=NULL) {
  if (!file.exists(outdir)) {
    dir.create(outdir, recursive=TRUE)
  }
  
  DATE <- format(Sys.time(), "%m%d%Y")
  
  if (is.null(baseName)) {
    fn <- paste0(test, "_", DATE, ".csv")
    qfn <- paste0(test, "_", DATE, ".png")
  } else {
    fn <- paste0(baseName, "_", test, "_", DATE, ".csv")
    qfn <- paste0(baseName, "_", test, "_", DATE, ".png")
  }
  
  # write the csv file version of the results
  rf <- file.path(outdir, fn)
  o <- order(results$p)
  write.table(results[o,], file=rf, sep=',', col.names=TRUE, row.names=FALSE, quote=TRUE)
  
  # write the qq plot
  qf  <- file.path(outdir, qfn)
  png(qf)
  qqunif(results$p, main=test, ci=TRUE)
  dev.off()
  
  # TBD: manhattan plot of results by chromosome
  return(invisible(NULL))
}


createTarballs <- function(files, ntar=1L, tarfile) {
  FILES <- na.omit(files)
  n <- length(FILES)
  
  outdir <- dirname(FILES[1])
  
  tardir <- dirname(tarfile)  
  if (!file.exists(tardir)) {
    dir.create(tardir, recursive=TRUE)
  }
  
  if (ntar == 1L) {
    tf <- tarfile
    f2tar <- paste(basename(FILES[1:n]), collapse=" ")
    tar_cmd <- paste("tar -cf", tf, "-C", outdir, f2tar, sep=" ")
    system(tar_cmd)  
  } else {
    a <- floor(n/ntar)
    for (i in 1:ntar) {
      PART <- paste("_part", i, ".tar", sep='')
      tf <- sub("\\.tar$", PART, tarfile)
      f2tar <- paste(basename(FILES[(((i-1)*a)+1):(i*a)]), collapse=" ")
      tar_cmd <- paste("tar -cf", tf, "-C", outdir, f2tar, sep=" ")
      system(tar_cmd)
    }
  }
}


saveFile <- function(obj, outdir, bn, suffix, type) {
  if (type == "RData") {
    filename <- paste(paste(bn, suffix, sep='_'), "RData", sep=".")
    message("...Saving ", suffix, " Object: ", filename)
    out <- file.path(outdir, filename)
    save(obj, file=out, compress="bzip2")    
  } else if (type == "csv") {
    filename <- paste(paste(bn, suffix, sep='_'), "csv", sep=".")
    message("...Writing ", suffix, " Object: ", filename)
    out <- file.path(outdir, filename)
    write.csv(obj, file=out)    
  }
}


Robj2txt <- function(x, file = paste(deparse(substitute(x)), ".txt", sep = "")) {
  tmp.wid = getOption("width")  # save current width
  options(width = 10000)        # increase output width
  sink(file)                    # redirect output to file
  print(x)                      # print the object
  sink()                        # cancel redirection
  options(width = tmp.wid)      # restore linewidth
  return(invisible(NULL))       # return (nothing) from function
} 

