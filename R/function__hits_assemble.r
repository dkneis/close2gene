
# for a single hits
assemble_one <- function(
  qstart, qend, qlen, qseq, sstart, send, slen, sseq, sstrand, gaps
) {
  # validation of inputs
  stopifnot(nchar(sseq) == slen)
  stopifnot(min(sstart, send) >= 1)
  stopifnot(max(sstart, send) <= slen)
  stopifnot(sstrand %in% c("plus","minus"))
  stopifnot(nchar(qseq) == qlen)
  stopifnot(qstart >= 1)
  stopifnot(qend <= qlen)
  stopifnot(qstart < qend)  # query sequence is always ordered this way
  # deal with hits on complementary strand
  if (sstrand == "minus") {
    sseq <- seq_reverseComplement(sseq)
    sstart.new <- slen - send + 1
    send.new <- slen - sstart + 1
    sstart <- send.new  # is correct
    send <- sstart.new  # is correct
  }
  # init output
  aseq <- NA
  amsg <- ""
  if ((gaps > 0) || ((send - sstart) != (qend - qstart))) {
    amsg <- "cannot assemble this because of gaps in alignment"
  } else {
    if ((qstart == 1) && (qend < qlen)) {
      if ((sstart == 1) && (send < slen)) {
        amsg <- paste("cannot assemble this case because overhang of",
          "both query and subject appear to the right of the aligned part")
      } else if ((sstart > 1) && (send == slen)) {
        aseq <- paste0(
          substr(sseq, 1, sstart - 1),
          substr(sseq, sstart, send),
          substr(qseq, qend + 1, qlen)
        )
        stopifnot(nchar(aseq) == qlen + slen - (send-sstart+1))
        if (sstrand == "minus") {
          aseq <- seq_reverseComplement(aseq)
        }
      } else {
        amsg <- paste("bad range of indices for subject sequence:",
          "aligned part must cover one (and only one) end")
      }
    } else if ((qstart > 1) && (qend == qlen)) {
      if ((sstart == 1) && (send < slen)) {
        aseq <- paste0(
          substr(qseq, 1, qstart - 1),
          substr(sseq, sstart, send),
          substr(sseq, send + 1, slen)
        )
        stopifnot(nchar(aseq) == qlen + slen - (send-sstart+1))
        if (sstrand == "minus") {
          aseq <- seq_reverseComplement(aseq)
        }
      } else if ((sstart > 1) && (send == slen)) {
        amsg <- paste("cannot assemble this case because overhang of",
          "both query and subject appear to the left of the aligned part")
      } else {
        amsg <- paste("bad range of indices for subject sequence:",
          "aligned part must cover one (and only one) end")
      }
    } else {
        amsg <- paste("bad range of indices for query sequence:",
          "aligned part must cover one (and only one) end")
    }
  }
  list(aseq=aseq, amsg=amsg)
}

# for all records in table of hits
assemble_all <- function(x) {
  # assemble
  x <- cbind(x, aseq=NA, amsg="")
  for (i in 1:nrow(x)) {
    out <- assemble_one(
      qstart=x[i,"qstart"], qend=x[i,"qend"], qlen=x[i,"qlen"], qseq=x[i,"qseq"],
      sstart=x[i,"sstart"], send=x[i,"send"], slen=x[i,"slen"], sseq=x[i,"sseq"],
      sstrand=x[i,"sstrand"], gaps=x[i,"gaps"]
    )
    x[i,"aseq"] <- out["aseq"]
    x[i,"amsg"] <- out["amsg"]
  }
  x <- x[!is.na(x[,"aseq"]),]
  if (nrow(x) > 0) {
    # remove shorter assemblies fully contained in longer ones
    x <- x[order(nchar(x[,"aseq"]), decreasing=TRUE),]
    out <- x[1,]
    if (nrow(x) > 1) {
      for (i in 2:nrow(x)) {
        if (! any(grepl(x[1:(i-1),"aseq"], pattern=x[i,"aseq"], fixed=TRUE))) {
          out <- rbind(out, x[i,])
        }     
      }
    }
    message(paste0("dropped ",nrow(x)-nrow(out)," redundant case(s)"))
    x <- out
  }
  x
}



