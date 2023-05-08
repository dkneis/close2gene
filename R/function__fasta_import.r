
fasta_import <- function(f, selection=NULL, blanks=FALSE) {
  tf <- tempfile()
  append <- FALSE
  con <- file(f, "r")
  while (TRUE) {
    hd <- readLines(con, 1)
    if (length(hd) == 0) {
      break  # done with file
    }
    if (substr(hd, 1, 1) != ">")
      stop("fasta file has line breaks or is corrupted")
    sq <- readLines(con, 1)
    hd <- substr(hd, 2, nchar(hd))
    if (!blanks && grepl(hd, pattern=" "))
      stop("whitespace in fasta headers not supported")
    if (is.null(selection) || (hd %in% selection)) {
      out <- write(paste0(hd,"\t",sq,"\n"), file=tf, append=append)
      append <- TRUE
    }
  }
  close(con)
  if (!append) {
    stop(paste0("no sequences imported from ",f))
  }
  out <- read.table(tf, sep="\t", header=FALSE)
  names(out) <- c("id","sequence")
  rm(tf)
  out
}
