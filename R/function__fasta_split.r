
fasta_split <- function(fasta.in, selection, fasta.out.selec, fasta.out.other) {
  stopifnot(length(selection) > 0)
  # add ">"
  selection <- paste0(">",selection)
  # walk through file
  con.in <- file(fasta.in, "r")
  con.out.selec <- file(fasta.out.selec, "w")
  con.out.other <- file(fasta.out.other, "w")
  count.selec <- 0
  count.other <- 0
  while (TRUE) {
    hd <- readLines(con.in, 1)
    if (length(hd) == 0) {
      break  # done with file
    }
    if (substr(hd, 1, 1) != ">")
      stop("fasta file has line breaks or is corrupted")
    sq <- readLines(con.in, 1)
    if (hd %in% selection) {
      write(x=hd, file=con.out.selec, append=(count.selec != 0))
      write(x=sq, file=con.out.selec, append=TRUE)
      count.selec <- count.selec + 1
    } else {
      write(x=hd, file=con.out.other, append=(count.other != 0))
      write(x=sq, file=con.out.other, append=TRUE)
      count.other <- count.other + 1
    }
  }
  close(con.in)
  close(con.out.selec)
  close(con.out.other)
  c(count.selec=count.selec, count.other=count.other)
}
