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
        help="The input fasta to align (required)."
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

  if (length(dna) <= 1) {
    writeXStringSet(dna, args$outfile)
    quit(save = "no", status = 0, runLast = FALSE)
  }

  # De-replicate the input sequences.
  u_dna <- unique(dna)
  index <- match(dna, u_dna)

  # Do alignment of the dereplicated sequences
  sink(stderr(), type = "output")
  U_DNA <- AlignSeqs(u_dna, processors=1)

  # Re-replicate the sequences.
  DNA <- U_DNA[index]
  names(DNA) <- names(dna)

  # Write to file.
  writeXStringSet(DNA, args$outfile)
}

main(args)
