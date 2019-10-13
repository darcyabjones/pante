#!/usr/bin/env Rscript

VERSION = "0.0.1"

suppressPackageStartupMessages(library("optparse"))
suppressWarnings(suppressPackageStartupMessages(library("Biostrings")))
suppressWarnings(suppressPackageStartupMessages(library("DECIPHER")))

option_list <- list(
    make_option(
        c("-i", "--infile"),
        type="character",
        action="store",
        help="The input aligned fasta (required)."
    ),
    make_option(
        c("-o", "--outfile"),
        type="character",
        action="store",
        help="The output file to write to (required)."
    ),
    make_option(
        "--version",
        type="logical",
        action="store_true",
        default=FALSE,
        help="Print version and exit.",
    )
)

parser <- OptionParser(
    usage = "%prog --infile in.fasta --outfile out.fasta",
    option_list = option_list
)

args <- parse_args(parser)

log_stderr <- function(...) {
  cat(sprintf(...), sep='', file=stderr())
}

quit_with_err <- function(...) {
  log_stderr(...)
  quit(save = "no", status = 1, runLast = FALSE)
}

validate_file <- function(path) {
  if (is.null(path)) {
    quit_with_err("Please provide required file")
  }
}


main <- function(args) {
  if (args$version) {
    cat(VERSION, file=stdout())
    quit(save = "no", status = 0, runLast = FALSE)
  }

  validate_file(args$infile)
  validate_file(args$outfile)

  dna <- readDNAStringSet(args$infile)

  mat <- consensusMatrix(dna, baseOnly=TRUE)[1:4,]
  mat <- rbind("N" = 0, mat)
  seq <- paste(
    apply(
      mat,
      MARGIN = 2,
      FUN=function(x) {names(which.max(x))}
    ),
    sep = "",
    collapse = ""
  )

  fam_name <- basename(tools::file_path_sans_ext(args$infile))
  names(seq) <- fam_name
  consensus <- DNAStringSet(seq)

  # Write to file.
  writeXStringSet(consensus, args$outfile)
}

main(args)
