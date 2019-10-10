#!/usr/bin/env python3

import sys
import argparse

from Bio import SeqIO
from Bio import AlignIO

from Bio.Alphabet import Gapped, SingleLetterAlphabet


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description=""" Extracts all sequences from a stockholm msa file into a
        single fasta without gaps.
        """
    )

    parser.add_argument(
        "infile",
        default=sys.stdin,
        type=argparse.FileType('r'),
        help="Input stockholm files. Use '-' for stdin.",
    )

    parser.add_argument(
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help="Output fasta file path. Default stdout.",
    )

    return parser.parse_args(args)


def main():
    args = cli(sys.argv[0], sys.argv[1:])
    seqs = []
    alignments = AlignIO.parse(
        args.infile,
        format="stockholm",
        alphabet=Gapped(SingleLetterAlphabet(), "-")
    )

    i = 1
    for alignment in alignments:
        for seq in alignment:
            seq.id = f"seq{i}"
            seq.name = f"seq{i}"
            seq.description = f"seq{i}"
            seq.seq = seq.seq.ungap()
            seqs.append(seq)
            i += 1

    SeqIO.write(seqs, args.outfile, "fasta")
    return


if __name__ == "__main__":
    main()
