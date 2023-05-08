
seq_unique <- function(x, colSeq) {
  # longest sequence first
  x <- x[order(nchar(x[,colSeq]), decreasing=TRUE),]
  # save longest
  keep <- x[1,]
  # decide on removal of others
  if (nrow(x) > 1) {
    for (i in 2:nrow(x)) {
      dup_plus <- any(grepl(keep[,colSeq], pattern=x[i,colSeq], fixed=TRUE))
      dup_minus <- any(grepl(keep[,colSeq], pattern=seq_reverseComplement(x[i,colSeq]), fixed=TRUE))
      if (!dup_plus && !dup_minus) {
        keep <- rbind(keep, x[i,])
      }
    }
  }
  keep
}
