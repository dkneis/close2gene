rm(list=ls())
options(stringsAsFactors=FALSE)

#################################################
# hard-coded global settings
#################################################

separ <- "+" # character to separate the IDs of individual
             # reads in the header of assembled contigs

#################################################
# import functions
#################################################

funcs <- list.files(path=".", pattern="^function__")
if (length(funcs) == 0) {
  stop("failed to load functions from current directory")
}
invisible(sapply(funcs, source))

#################################################
# get command line args as list
#################################################

arg <- commandArgs(trailingOnly=TRUE)
arg <- arg[arg != "--args"]
if (!all(grepl(arg, pattern="^[^=]+[=][^=]+$" ))) {
  stop("arguments after '--args' must take the form 'name=value'")
}
n <- gsub(arg, pattern="[=]([^=]+)$", replacement="")
v <- gsub(arg, pattern="^([^=]+)[=]", replacement="")
#print(n)
#print(v)
arg <- list(v)
arg <- lapply(v, type.convert, as.is=TRUE)
names(arg) <- n
message("running with the following user inputs:")
for (i in 1:length(arg)) {
message(paste0("  ",names(arg[i]),"=",arg[[i]]))
}
rm(n,v)

getArg <- function(name) {
  x <- arg[[name]]
  if (is.null(x)) {
    stop(paste0("missing assignment '",name,"=...' in command line"))
  }
  x
}

#################################################
# assign and check arguments
#################################################

seqInitPool <- getArg("pool")
if (!file.exists(seqInitPool))
  stop(paste0("file provided as 'pool' not found: ",seqInitPool))

seqInitSeed <- getArg("seed")
if (!file.exists(seqInitSeed))
  stop(paste0("file provided as 'seed' not found: ",seqInitSeed))

odir <- getArg("odir")
if (!dir.exists(odir))
  stop(paste0("output directory specified via 'odir' not found: ",odir))

prefix <- getArg("prefix")
if (grepl(prefix, "/"))
  stop("cannot use string specified as 'prefix' to name output files")

pool_compressed <- getArg("pool_compressed")
if (!is.logical(pool_compressed))
  stop("value of 'pool_compressed' must be TRUE or FALSE")

newlines <- getArg("newlines")
if (!is.logical(newlines))
  stop("value of 'newlines' must be TRUE or FALSE")

cleanup <- getArg("cleanup")
if (!is.logical(cleanup))
  stop("value of 'cleanup' must be TRUE or FALSE")

max_rank <- getArg("max_rank")
if (!is.numeric(max_rank))
  stop("value of 'max_rank' must be numeric")
if (max_rank < 1)
  stop("value of 'max_rank' not in reasonable range")

min_pident <- getArg("min_pident")
if (!is.numeric(min_pident))
  stop("value of 'min_pident' must be numeric")
if ((min_pident < 50) || (min_pident > 100))
  stop("value of 'min_pident' not in reasonable range")

min_length <- getArg("min_length")
if (!is.numeric(min_length))
  stop("value of 'min_length' must be numeric")
if (min_length < 8)
  stop("value of 'min_length' not in reasonable range")

max_iter <- getArg("max_iter")
if (!is.numeric(max_iter))
  stop("value of 'max_iter' must be numeric")
if (max_iter < 1)
  stop("value of 'max_iter' not in reasonable range")

#################################################
# check availability of dependencies
#################################################

dependencies <- c("tar", "blastn", "makeblastdb")
if (Sys.info()["sysname"] == "Linux") {
  ok <- TRUE
  for (dep in dependencies) {
    if (system2("which", dep, stdout=FALSE) != 0) {
      ok <- FALSE    
    }
  }
  if (!ok)
    stop(paste0("missing dependencies; need: '",
      paste(dependencies, collapse=", "),"'"))
} else {
  message("Not running on Linux - Check for dependencies skipped!") 
}

#################################################
# unzip / copy sequences to be scanned
#################################################

if (pool_compressed) {
  message("extracting pool archive")
  cmd <- "tar"
  arg <- paste0("-xzf ",seqInitPool," --directory=",odir)
  stat <- system2(cmd, arg, stdout=FALSE, stderr=FALSE)
  if (stat != 0)
    stop(paste0("extraction of input archive failed with exit code ",stat))
  seqInitPool <- paste0(odir,"/",gsub(basename(seqInitPool), pattern="[.]tar[.]gz$", replacement=""))
} else {
  # making a copy is important because we delete the old pool after each iteration;
  # without this copying, the original input file would get lost
  message("copying pool")
  file.copy(seqInitPool, odir)
  seqInitPool <- paste0(odir,"/",basename(seqInitPool))
}
if (!file.exists(seqInitPool))
  stop(paste0("could not find pool in expected path: ",seqInitPool))

#################################################
# relabel query sequences
#################################################

if (newlines) {
  message("removing newlines from pool")
  fasta_nobreaks(
    fasta.in=seqInitPool,
    fasta.out=paste0(seqInitPool,".nobreaks")
  )
  invisible(file.remove(seqInitPool))
  invisible(file.rename(paste0(seqInitPool,".nobreaks"), seqInitPool))
}

message("relabeling sequences in pool")
fasta_relabel(
  fasta.in=seqInitPool,
  fasta.out=paste0(seqInitPool,".relabeled")
)
invisible(file.remove(seqInitPool))
invisible(file.rename(paste0(seqInitPool,".relabeled"), seqInitPool))

#################################################
# copy seed file
#################################################

message("copying seed")
fname <- paste0(odir,"/",prefix,".iter",1,"_seed.fasta")
invisible(file.copy(seqInitSeed, fname))
seqInitSeed <- fname

if (newlines) {
  message("removing newlines from seed")
  fasta_nobreaks(
    fasta.in=seqInitSeed,
    fasta.out=paste0(seqInitSeed,".nobreaks")
  )
  invisible(file.remove(seqInitSeed))
  invisible(file.rename(paste0(seqInitSeed,".nobreaks"), seqInitSeed))
}

# define outputs to be generated by blast using tabular output (outfmt 6)
# see "https://www.metagenomics.wiki/tools/blast/blastn-output-format-6" for defaults
# note that we add "qlen" and "slen" to also get the total length of the involved sequences;
# we need as least 'slen' for the analysis
blastOutCols <- c("qaccver","saccver","pident","length","mismatch",
  "gaps","qlen","qstart","qend","slen","sstart","send","sstrand",
  "evalue","bitscore")

#################################################
# start iteration
#################################################

for (iter in 1:max_iter) {
  
  message(paste0("iteration ",iter))
  
  if (iter == 1) {
    seqPool <- seqInitPool
    seqSeed <- seqInitSeed
  } else {
    seqPool <- seqNewPool
    seqSeed <- seqNewSeed
  }

  #################################################
  # note: We create the database to feed the subjects into
  #       blastn via the -db argument instead of the -subject
  #       argument. Memory consumption is reduced dramatically.
  message("creating database from seeds")
  dbSeed <- paste0(odir,"/",prefix,"_seedDB")
  cmd <- "makeblastdb"
  arg <- paste0(
    " -dbtype nucl -input_type fasta -hash_index ",
    " -in ",seqSeed,
    " -out ",dbSeed
  )
  stat <- system2(cmd, arg, stdout=FALSE, stderr=FALSE)
  if (stat != 0)
    stop(paste0("database creation failed with exit code ",stat))

  #################################################
  message("scanning for hits using blast")
  tblHits <- paste0(odir,"/",prefix,".iter",iter,"_hitsRaw.txt")
  cmd <- "blastn"
  arg <- paste0(
    " -query ",seqPool,
    " -db ",dbSeed,
    " -use_index false",
    " -outfmt \"6 ",paste(blastOutCols, collapse=" "),"\"",
    " -out ",tblHits
  )
  stat <- system2(cmd, arg, stdout="", stderr="")
  if (stat != 0)
    stop(paste0("blast search failed with exit code ",stat))
  
  #################################################
  message("reading hits")
  has.hits <- TRUE
  h <- tryCatch(read.table(tblHits, header=FALSE, sep="\t"),
    error = function(e) { has.hits <<- FALSE}
  )
  if (has.hits) {
    names(h) <- blastOutCols
  }

  if (has.hits) {
    message("filtering hits")
    h <- hits_filter(
      x= h,
      initial= (iter == 1),
      max_rank= max_rank,
      min_pident= min_pident,
      min_length= min_length
    )
    if (nrow(h) == 0) {
      message("no hits passed the filter")
      has.hits <- FALSE
    } else {
      write.table(x=h, file=paste0(odir,"/",prefix,".iter",iter,"_hitsFiltered.txt"), sep="\t",
        col.names=TRUE, row.names=FALSE, quote=FALSE)
    }
  }

  if (has.hits) {
    message("getting affected seed sequences and cleaning up")
    tmp <- fasta_import(f=seqSeed, selection=unique(h[,"saccver"]))
    h <- merge(x=h, y=tmp, by.x="saccver", by.y="id")
    names(h)[names(h) == "sequence"] <- "sseq"
    rm(tmp)
    #file.remove(seqSeed)

    message("getting and removing affected sequences from pool")
    seqNewPool <- paste0(odir,"/",prefix,".iter",iter+1,"_pool.fasta")
    tf <- tempfile()
    stat <- fasta_split(
      fasta.in=seqPool,
      selection=unique(h[,"qaccver"]),
      fasta.out.selec=tf,
      fasta.out.other=seqNewPool
    )
    if (any(stat == 0))
      stop("problem with splitting")
    tmp <- fasta_import(tf, selection=NULL)
    h <- merge(x=h, y=tmp, by.x="qaccver", by.y="id")
    names(h)[names(h) == "sequence"] <- "qseq"
    rm(tmp)
    invisible(file.remove(tf))
  }

  has.assembled <- TRUE
  if (has.hits) {
    message("performing assembly")
    h <- assemble_all(h)
    message(paste("successful assemblies: ",nrow(h)))
    if (nrow(h) == 0) {
      has.assembled <- FALSE
    } else {
      message("generating new seed")
      seqNewSeed <- paste0(odir,"/",prefix,".iter",iter+1,"_seed.fasta")
      h[,"saccver"] <- paste(h[,"saccver"], h[,"qaccver"], sep=separ)
      h <- seq_unique(h, colSeq="aseq")  # duplicate removal
      fasta_export(f=seqNewSeed, x=h, col_id="saccver", col_seq="aseq")
      rm(h)
    }
  }

  message("cleaning up temporary files")
  file.remove(seqPool)
  f <- paste0(dbSeed, c(".nhd",".nhi",".nhr",".nin",".nog",".nsd",".nsi",".nsq"))
  sapply(f, file.remove)

  if (!has.hits) {
    message(paste0("no hits found in iteration ",iter," of ",max_iter))
    break
  }

  if (!has.assembled) {
    message(paste0("nothing assembled in iteration ",iter," of ",max_iter))
    break
  }

}
if (has.hits && file.exists(seqNewPool)) {
  invisible(file.remove(seqNewPool))
}

message("creating summary files")
summarizeAssemblies(
  dir=odir,
  fasta.full=paste0(odir,"/",prefix,".summary_full.fasta"),
  fasta.ext=paste0(odir,"/",prefix,".summary_ext.fasta"),
  pre=paste0(prefix,".iter"),
  suf="_seed.fasta",
  separ="+"
)

if (cleanup) {
  message("removing intermediate outputs")
  f <- list.files(path=odir, pattern=paste0("^",prefix,"[.]iter[0-9]+_"), full.names=TRUE)
  invisible(file.remove(f))
}

message("done")
