#!/usr/bin/env python3

import sys
import argparse
import os

from Bio import AlignIO

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
        alignments = AlignIO.read(infile, format="fasta")

        id_ = os.path.split(os.path.splitext(infile.name)[0])[-1]

        args.outfile.write("#=GF ID {}\n".format(id_))
        AlignIO.write(alignments, args.outfile, "stockholm")

    return


if __name__ == "__main__":
    main()

