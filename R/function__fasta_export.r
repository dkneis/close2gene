
fasta_export <- function(f, x, col_id, col_seq) {
  stopifnot(all(c(col_id,col_seq) %in% names(x)))
  append <- FALSE
  con <- file(f, "w")
  for (i in 1:nrow(x)) {
    write(file=con, x=paste0(">",x[i,col_id],"\n",
      x[i,col_seq]), append=append)
    append <- TRUE
  }
  close(con)
}
