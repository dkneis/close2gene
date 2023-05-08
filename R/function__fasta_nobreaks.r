
fasta_nobreaks <- function(fasta.in, fasta.out) {
  con.in <- file(fasta.in, "r")
  con.out <- file(fasta.out, "w")
  append <- FALSE
  sq <- ""
  while (TRUE) {
    x <- readLines(con.in, 1)
    if (length(x) == 0) {
      if (sq != "") {
        write(x=sq, file=con.out, append=TRUE)
      }
      break  # done with file
    }
    if (substr(x, 1, 1) == ">") {
      if (sq != "") {
        write(x=sq, file=con.out, append=TRUE)
      }
      write(x=x, file=con.out, append=append)
      append <- TRUE
      sq <- ""
    } else {
      sq <- paste0(sq, x)
    }
  }
  close(con.in)
  close(con.out)
}
