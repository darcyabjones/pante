#!/usr/bin/env Rscript

VERSION = "0.0.1"

suppressPackageStartupMessages(library("optparse"))
suppressWarnings(suppressPackageStartupMessages(library("DECIPHER")))
suppressWarnings(suppressPackageStartupMessages(library("Biostrings")))
data(PFASUM)

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
        default="-",
        help="The output file to write to (required)."
    ),
    make_option(
        c("-t", "--tmpdir"),
        type="character",
        action="store",
        default=tempdir(),
        help="Where to store temporary files if using stdin or stdout"
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

# From https://github.com/Bioconductor/Biostrings/issues/10
from_string_to_XString <- function(x, Class, alphabet) {
    stopifnot(isSingleString(x),
              isSingleString(Class),
              is.character(alphabet), !anyNA(alphabet))
    keys <- vals <- as.integer(charToRaw(paste0(alphabet, collapse="")))
    lkup <- Biostrings:::buildLookupTable(keys, vals)
    .Call2(
      "new_XString_from_CHARACTER",
      Class, x, 1L, nchar(x), lkup,
      PACKAGE="Biostrings"
    )
}


main <- function(args) {

  if (args$version) {
    cat(VERSION, file=stdout())
    quit(save = "no", status = 0, runLast = FALSE)
  }

  if (args$infile == "-") {
    infilename = tempfile(tmpdir = args$tmpdir, fileext = ".fasta")
    sin <- readLines("stdin")
    writeLines(sin, con=infilename)
  } else {
    infilename = args$infile
  }

  if (args$outfile == "-") {
    outfilename = tempfile(tmpdir = args$tmpdir, fileext = ".fasta")
  } else {
    outfilename = args$outfile
  }

  validate_file(infilename)
  validate_file(outfilename)

  seqs <- readAAStringSet(infilename)

  seqs <- endoapply(
    seqs, 
    function(x) {
      from_string_to_XString(
        x,
        "AAString",
        AA_ALPHABET
      )
    }
  )

  if (args$infile == "-" && file.exists(infilename)) {
    dummy <- file.remove(infilename)
  }

  # De-replicate the input sequences.
  u_seqs <- unique(seqs)
  index <- match(seqs, u_seqs)

  # No need to align identical or single sequences.
  if (length(u_seqs) <= 1) {
    writeXStringSet(seqs, args$outfile)
    quit(save = "no", status = 0, runLast = FALSE)
  }

  # Align dereplicated sequences
  tryCatch(
    expr = {
      aligned <- AlignSeqs(
        u_seqs,
        iterations=2,
        refinements=2,
        substitutionMatrix=PFASUM[,,50],
        verbose=FALSE,
        processors=1
      )
    },
    error = function(e) {
      warning(paste(names(seqs), collapse="\n"))
      stop(e)
    }
  )

  # Re-replicate the sequences.
  SEQS <- aligned[index]
  names(SEQS) <- names(seqs)

  writeXStringSet(SEQS, outfilename)

  if (args$outfile == "-") {
    sout <- readLines(outfilename)
    writeLines(sout)

    if (file.exists(outfilename)) {
      dummy <- file.remove(outfilename)
    }
  }


}

main(args)
