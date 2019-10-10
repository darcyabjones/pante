#!/usr/bin/env python3

import sys
import argparse
import os

from Bio import AlignIO

from Bio.Alphabet import Gapped, SingleLetterAlphabet


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description=""" Extracts all sequences from a fasta msa file into a
        stockholm alignment file.
        """
    )

    parser.add_argument(
        "infiles",
        type=argparse.FileType('r'),
        nargs="+",
        help="Input fasta files.",
    )

    parser.add_argument(
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help="Output stockholm file. Default stdout.",
    )

    return parser.parse_args(args)


def main():
    args = cli(sys.argv[0], sys.argv[1:])

    for infile in args.infiles:
        alignments = AlignIO.read(
            infile,
            format="fasta",
            alphabet=Gapped(SingleLetterAlphabet(), "-")
        )

        id_ = os.path.split(os.path.splitext(infile.name)[0])[-1]

        fmt = alignments.format("stockholm").split("\n", maxsplit=1)
        args.outfile.write(fmt[0])
        args.outfile.write(f"\n#=GF ID {id_}\n")
        args.outfile.write(fmt[1])
    return


if __name__ == "__main__":
    main()
