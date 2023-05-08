rm(list=ls())

source("../../../R/function__seq_reverseComplement.r")

ifile <- "pB10.fasta"  # must be single-fasta without newlines in the sequence!
readLength <- 260
meanCoverage <- 5
ofile <- paste0("../pB10_len",readLength,"_cov",meanCoverage,".fasta")

con <- file(ifile, "r")
header <- readLines(con=con, n=1)
genome <- readLines(con=con, n=1)
close(con)

numReads <- nchar(genome) / readLength * meanCoverage
starts <- runif(numReads, min=1, max=nchar(genome)-readLength)
fn <- function(i) {
  write(file=con, x=paste0(">seq_",i), ncolumns=1)
  seq <- substr(genome, starts[i], starts[i]+readLength)
  if (runif(1) > 0.5)
    seq <- seq_reverseComplement(seq)
  write(file=con, x=seq, ncolumns=1)
}
con <- file(ofile, "w")
sapply(1:length(starts), fn)
close(con)
