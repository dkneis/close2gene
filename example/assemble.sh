#!/bin/bash

cd ../R
Rscript --vanilla main.r --args\
  pool=../example/in/pB10_len260_cov5.fasta\
  seed=../example/in/sul1.fasta\
  odir=../example/out\
  prefix=sul1_in_pB10\
  pool_compressed=FALSE\
  newlines=FALSE\
  cleanup=TRUE\
  max_rank=3\
  min_pident=95\
  min_length=50\
  max_iter=3
