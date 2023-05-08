
seq_reverseComplement <- function(x) {
  # reverse
  x <- paste0(rev(unlist(strsplit(x, split=""))),collapse="")
  # complement
  x <- toupper(x)
  x <- gsub(x, pattern="A", replacement="t", fixed=TRUE)
  x <- gsub(x, pattern="T", replacement="a", fixed=TRUE)
  x <- gsub(x, pattern="C", replacement="g", fixed=TRUE)
  x <- gsub(x, pattern="G", replacement="c", fixed=TRUE)
  x <- toupper(x)
  x
}
