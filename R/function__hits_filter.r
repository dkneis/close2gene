
hits_filter <- function(x, initial, max_rank, min_pident, min_length) {
  stopifnot(c("qaccver","saccver","evalue","pident","length") %in% names(x))
  message(paste0("records in input table: ",nrow(x)))
  if (initial) {
    # in the very first pass, we want to find all sequences in the pool that match
    # a seed sequence but each sequence from the poll should only be linked with the
    # best-matching seed
    x <- by(x, x[,"qaccver"], function(z) {z[which.min(z[,"evalue"]),]})
    x <- do.call(rbind, x)
    message(paste0("records remaining after selection of best-matching seed: ",nrow(x)))
  } else {
    # in the phase of extension, we want each (grown) seed to be extended by the best-matching
    # sequence(s) from the pool; here, we can tolerate multiple matches (leading to alternative assembly paths)
    x <- by(x, x[,"saccver"], function(z) {z[which(rank(z[,"evalue"]) <= max_rank),]})
    x <- do.call(rbind, x)
    message(paste0("records remaining after selection of best-matching extension: ",nrow(x)))
  }
  x <- x[with(x, (pident >= min_pident) & (length >= min_length)),]
  message(paste0("records remaining after quality filter: ",nrow(x)))
  x
}

