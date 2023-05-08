
fasta_relabel <- function(fasta.in, fasta.out) {
  con.in <- file(fasta.in, "r")
  con.out <- file(fasta.out, "w")
  append <- FALSE
  i <- 0
  while (TRUE) {
    i <- i + 1
    hd <- readLines(con.in, 1)
    if (length(hd) == 0) {
      break  # done with file
    }
    if (substr(hd, 1, 1) != ">")
      stop("fasta file has line breaks or is corrupted")
    hd <- paste0(">",i)
    sq <- readLines(con.in, 1)
    write(x=hd, file=con.out, append=append)
    append <- TRUE
    write(x=sq, file=con.out, append=TRUE)
  }
  close(con.in)
  close(con.out)
}
