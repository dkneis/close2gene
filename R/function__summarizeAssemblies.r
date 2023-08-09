
summarizeAssemblies <- function(
  dir,
  fasta.full,
  fasta.ext,
  pre,  # file part preceding the iteration number
  suf,  # file part after the iteration numner
  separ="+", # character separating the read IDs in the header of the assembly
  silent=FALSE
) {
  stopifnot(dir.exists(dir))
  pre <- gsub(pre, pattern=".", replacement="[.]", fixed=TRUE)
  suf <- gsub(suf, pattern=".", replacement="[.]", fixed=TRUE)
  # get list of files
  files <- list.files(path=dir, pattern=paste0(pre,"[0-9]+",suf), full.names=TRUE)
  if (length(files) > 0) {
    # collect from all iterations
    x <- NULL
    for (f in files) {
      iter <- as.integer(gsub(basename(f), pattern=paste0(pre,"([0-9]+)",suf), replacement="\\1"))
      x <- rbind(x, cbind(
        iteration=iter, fasta_import(f, selection=NULL)
      ))
    }
    if (!identical(sort(unique(x[,"iteration"])), 1:length(files))) {
      stop(paste0("collection of files in directory '",dir,"' is apparently incomplete"))
    }
    # drop shorter assemblies being contained in longer ones
    x <- x[order(nchar(x[,"sequence"]), decreasing=TRUE),]
    out <- x[1,]
    if (nrow(x) > 1) {
      for (i in 2:nrow(x)) {
        if (! any(grepl(x[1:(i-1),"sequence"], pattern=x[i,"sequence"], fixed=TRUE))) {
          out <- rbind(out, x[i,])
        }     
      }
    }
    if (!silent) {
      message(paste0("number of result sequences reduced from ",nrow(x)," to ",nrow(out)))
    }
    # export full sequences containing the seed
    fasta_export(f=fasta.full, x=out, col_id="id", col_seq="sequence")
    # export only the sequence which has beed appended to the seed
    seeds <- x[x[,"iteration"] == 1, c("id","sequence")]
    
    # this function is needed since grepl / gsub don't work properly
    # if the pattern argument has several thousand characters
    delete_seed_at_margin <- function(ass, seed) {
      if (nchar(seed) <= nchar(ass)) {
        if (substr(ass, 1, nchar(seed)) == seed) {
          out <- paste0("N",substr(ass, nchar(seed)+1, nchar(ass)))
        } else if (substr(ass, nchar(ass) - nchar(seed) + 1, nchar(ass)) == seed) {
          out <- paste0(substr(ass, 1, nchar(ass) - nchar(seed)),"N")
        } else {
          stop("could not identify seed in assembled sequence")
        }
      } else {
        stop("could not identify seed in assembled sequence")        
      }
      out
    }
    
    for (i in 1:nrow(out)) {
      s <- gsub(out[i,"id"], pattern=paste0("[",separ,"].*$"), replacement="")
      k <- which(seeds[,"id"] == s)
      if (length(k) != 1) {
        stop("failed to determine seed ID from header of assembly")
      }
#      if (grepl(out[i,"sequence"], pattern=paste0("^",seeds[k,"sequence"]))) {
#        out[i,"sequence"] <- gsub(out[i,"sequence"], pattern=paste0("^",seeds[k,"sequence"]), replacement="N")
#      } else if (grepl(out[i,"sequence"], pattern=paste0(seeds[k,"sequence"],"$"))) {
#        out[i,"sequence"] <- gsub(out[i,"sequence"], pattern=paste0(seeds[k,"sequence"],"$"), replacement="N")
#      } else {
#        stop("failed to delete seed from assembled sequence")
#      }
      out[i,"sequence"] <- delete_seed_at_margin(out[i,"sequence"], seeds[k,"sequence"])
    }
    fasta_export(f=fasta.ext, x=out, col_id="id", col_seq="sequence")
  } else {
    stop("no input files with matching pattern")
  }
}

#source("function__fasta_import.r")
#source("function__fasta_export.r")
#summarizeAssemblies(
#  dir="/home/dkneis/tudd/dev/bioinfo/seedcontext/test/K12/out",
#  fasta.full="summary_full.fasta",
#  fasta.ext="summary_ext.fasta"
#)
