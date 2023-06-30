# About this software

## What is does
This repository contains R source code to perform a **seed-based assembly
of DNA** reads to identify the flanking region of a particular sequence
of interest (e.g. the context of an antibiotic resistance gene). It was
originally developed to identify the context of *dfrB* genes in
shotgun metagenomic data of environmental origin by [Kneis et al., 2023, ISMEJ](10.1038/s41396-023-01460-7).

The major practical advantage of seed-based assembly is the **very low RAM
consumption**. Classical, graph-based assemblers may easily require hundreds
of GB (or even several TB) of memory in the processing of environmental
metagenomes. With the seed-based approach, RAM usage is typically not a concern
and several input files may be processed in parallel with standard resources of,
e.g., 20 GB of RAM and 8 CPUs.

Note that **similar software** was developed in the past by other groups. You
might want to consider, e.g., the software described in the 2008 paper of
Sobreira and Gruber
[DOI: 10.1093/bioinformatics/btn283](https://doi.org/10.1093/bioinformatics/btn283).

## Requirements

The source is written in the [R language](https://www.r-project.org/) and,
consequently, you need a working installation of R to run it.

The R code contains calls to a few external programs, all of which are available
free of charge. Those are:

- `blastn` : the well known tool performing alignment of nucleotide sequences
  which can be obtained from [this NCBI website](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#blast-executables).

- `makeblastdb` : a helper program to pre-process input for `blastn`. If you
  installed `blastn`, you probably have this tool on board already.

- `tar` : an archiving tool used to compress / unzip files which can be obtained,
  for example, [from this GNU website](https://www.gnu.org/software/tar/). Users
  of Linux probably have the tool installed.

Make sure that these programs are actually findable, i.e. their respective
installation directories are registered in the PATH environment variable.

So far, the software was exclusively tested on **Linux** platforms, namely Ubuntu
20.04 LTS. Other platforms may work but the developer does not intent to
provide specific support for any other platform. Why would you use something
else than Linux for serious computing in an academic context?

## License
The software is distributed under the Gnu General Public License. It thus
comes without any warranty - see the attached LICENSE file. This licence
applies to all source files in the repository.

## Contact
In case of questions, please contact the author, David Kneis, a member of the
TU Dresden-based Institute of Hydrobiology (firstname.lastname@tu-dresden.de).

# Contents of the repository

The source consists of a collection of R files organized as follows:

- `main.r` : This is the main program.

- `function__*.r` : Each of these files contains the source code of a specific
  function called by the main program. The contents of all files is automatically
  imported upon execution of `main.r`.

At a later stage of development, the source code might be re-organized in the
form of a classical R package.

# Usage

## Basic execution strategy

In contrast to typical R applications, the program is not intended to be run in
an interactive mode. You rather want to trigger execution from the command line
such that sequence analysis can be performed on a powerful server without a
graphical interface. The usual way of doing this is via the `Rscript` program
that comes with R.

Generally, execution could be invoked like below where `ARGS` is a placeholder
for any command line arguments that you need to provide to the program (see below
for details). Check the help of `Rscript` for the meaning of `--vanilla` and `--args`.

```
Rscript --vanilla main.r --args ARGS
```

Computations can be time consuming. If you work on a server, you may thus want
log out while computations are running or a disruption of the network might
trigger disconnection. Therefore, on a server, you most likely want to run the
program like so

```
nohup Rscript --vanilla main.r --args ARGS >LOGFILE 2>&1 &
```

to prevent hangup, redirect all error and standard outputs to LOGFILE, and
run the process in the background in order not to block the current terminal
until things are finished.

## Command line arguments

### Example call

Command line arguments must be passed in `key=value` format and a complete call
may look like so:

```
nohup Rscript --vanilla main.r --args pool=../metagenomes/SRR9999.fasta.tar.gz
  seed=genesOfInterest.fasta odir=out prefix=SRR9999 pool_compressed=TRUE
  newlines=FALSE cleanup=TRUE max_rank=3 min_pident=95 min_length=50
  max_iter=10 >log/SRR9999.log 2>&1 &
```

### Individual arguments

- `pool` : File path pointing to file in fasta format holding sequences from
  metagenomic (or possibly genomic) sequences. In the context of paired-end
  sequencing, the file contains the merged and quality-trimmed forward and
  reverse reads. Optionally, the file can be provided in compressed format
  (.tar.gz).
  
- `seed` : File path pointing to a file in fasta format containing the seed
  sequences whose flanking regions shall be assembled using the sequences in
  the `pool` file. Usually, this is a small file with just a few entries. It
  must be a plain fasta file without compression.

- `pool_compressed` : Must be either TRUE or FALSE (capitals required). If TRUE,
  it is assumed that the `pool` file is compressed (.tar.gz) and needs to be
  unzipped prior to processing by `blastn`.

- `newlines` : Must be either TRUE or FALSE (capitals required). If TRUE,
  newline characters can be present in input fasta files. If FALSE, the inputs
  are assumed to be without line breaks which speeds up processing. Use the
  `head` command to check with case applies to you.

- `odir` : The name of an existing directory where any outputs created by the
  program should appear. It is recommended that the directory is initially empty
  to prevent overwriting of possibly existing files.

- `out prefix` : A prefix to be used as the first part of output file names.
  Typically, you want to use (part of) the basename of the file used as `pool`.

- `cleanup` : Must be either TRUE or FALSE (capitals required). If TRUE, all
  intermediate files will be automatically deleted to not clutter the output
  directory. It may be helpful to set this to FALSE if you need to understand
  unexpected outputs or find the reason for bad performance.

- `max_rank` : Given a seed sequence (or an already built contig), the pool of
  (remaining) sequences may contain several candidate reads all of which
  show a good alignment with either end of the seed/contig. The program
  automatically sorts the candidates by the quality of alignment. With this integer
  argument, you can control how many of the candidates shall be considered 
  in subsequent assembly steps. If you set this to a number greater than 1, the
  assembler will produce multiple possible trajectories per seed sequence or 
  contig which results in a tree-like structure of results. Most likely, you
  want a number between 1 (high coverage WGS data) and something between 3 and
  5 (metagenomes). Note that the amount of output may increase exponentially
  as higher numbers are chosen!

- `min_pident` : Controls the selection of candidate reads for contig extension.
  A read from the (remaining) pool of sequences will only be considered if the
  percentage identify of the alignment (as found by `blastn`) is at least as
  high as this value. Reasonable values are most likely > 90 but this
  depends on the requested alignment length (see below).

- `min_length` : Works together with the argument `min_pident`. The argument lets
  you specify the minimum number of base pairs which an alignment must have such
  that a candidate read is considered as a reasonable contig extension.
  Specifically, only alignments at the ends of sequences are considered, so this
  is effectively a minimum length of sequence overlap. If you are tolerant with
  regard to `min_pident`, you probably want to be more strict here and request a
  longer value (or vice versa).

- `max_iter` : Contigs are extended in an iterative manner. In each iteration,
  only one read from the (remaining) pool of sequences is added at the end of
  every contig. The algorithm (and thus the program) naturally stops when no
  further extension of any contig is possible. Depending on your data, this may
  result in a very large number of iterations and, very likely, unacceptably long
  computation times. Therefore, it is recommended to let the algorithm stop
  after a certain number of iterations. Example: If your input reads are 250 bp
  long and you request a minimum overlap of 50 bp, you will, at the very best, see
  an extension of the contig by 200 bp in each iteration. If you were happy with
  contigh length of 2000 bp, a reasonable approximate choice for `max_iter`
  would be 10.

### Additional details on input files

The fasta files provided as inputs must currently adhere to the following
basic conventions:

- The fasta headers, i.e. the ID of individual sequences, must not contain
  whitespace as this cause trouble with the processing of `blastn` outputs.
  Any whitespace should be replace by an alternative character like, e.g.,
  an underscore.
  
- All fasta input files should use a consistent case, i.e. nucleotides should
  either be encoded as "ATCG" or "atcg". Do not use different cases in distinct
  input files.

## Parallel processing

The program currently does not use parallel processing internally, namely
`blastn` is called without the `-num_threads` argument.

If you want to run the assembly algorithm for multiple inputs in parallel, this
can still be organized with little effort. The simplest way is to put the names of all
your input files into the column of a spreadsheet. Then, using a spreadsheet
function like `CONCAT`, you create a derived column holding the complete
call, e.g. `Rscript --vanilla main.r --args ARGS`. Say you want to use 4 CPUs.
In that case, you copy 1/4 of the lines from the spreadsheet into 4 separate
shell script files and run those with a sequence of commands like:

```
nohup ./part1.sh >part1.log 2>&1 &
nohup ./part2.sh >part2.log 2>&1 &
nohup ./part3.sh >part3.log 2>&1 &
nohup ./part4.sh >part4.log 2>&1 &
```

You might want to write down the respective ID of the parent process displayed
by `nohup`. Without that ID it may be difficult to stop one or all of the
instances, if needed, because of the various child processes created.

# Minimum working example

## Objective

We consider the pB10 plasmid described by
[Schlueter et al., (2003)](https://doi.org/10.1099/mic.0.26570-0) whose complete
nucleotide sequence is deposited under [GenBank accession number AJ564903.1](https://www.ncbi.nlm.nih.gov/nuccore/AJ564903.1).
The plasmid hosts several antibiotic resistance genes, among which is *sul1*.
In this example, we want to identify the genetic context of the sulfonamide
resistance gene *sul1*.

All material related to this example is contained in the "example" folder.

## Input data

Assume that we have the sequenced plasmid DNA available as reads of 260 bp
length. This would be the typical situation after the merging of paired-end
reads generated by 2 x 150 bp sequencing.

For simplicity, in this example, we generated fake reads from the known complete
sequence of the plasmid by randomly picking subsequences of 260 bp length. We
perform the sampling such that the statistical coverage is 5, i.e.  each part of
the complete sequence is, on average, covered by 5 reads. Note that the generated
fake reads are, unlike true sequencing results, free of errors.

The fasta file containing these fake reads is what will be passed to the
assembler via the `pool` argument. The file is named "pB10_len260_cov5.fasta"
and the script for creating this file can be found in the folder "make_fake_reads".

For the `seed` argument, we create another fasta file containing the sequence
of the *sul1* gene. Thus, in this example, we consider a single seed only. The
respective file name is "sul1.fasta".

## Program call

The assembly is performed by executing the bash script "assemble.sh" which
contains the following:

```
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
```

Thus, we are providing uncompressed files, without newlines. Further, we
tell the assembler to perform 3 iterations and request a minimum of 95%
nucleotide identity together with a minimum alignment length of 50 bp in the
selection of candidate reads for contig extension.

## Output

After a few seconds, two outpul files should appear in the "out" subfolder.

- `sul1_in_pB10.summary_full.fasta` : This file contains the assembled
  contigs which include the seed sequences.
  
- `sul1_in_pB10.summary_ext.fasta` : This file contains the same contigs but
  the seed sequences were removed. Specifically, the seed sequences were replaced
  by a single "N" character. This makes it possible to see at which end of the
  contig the seed sequence was located.
  
The second output file, which lacks the seed sequences, is particularly suited
for further analysis of the flanking region of the seed sequence. It can be,
for example, feed into the [`blastx` program](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
to possibly identify the neighborhood of the seed.

For the example data set, a `blastx` search suggests that *sul1* is flanked by genes
coding for an acetyltransferase and a multi-drug efflux pump of the QacE family.
This is in perfect agreement with the (known) order of genes in the complete
pB10 plasmid. Thus, the assembly was obviously successful.
