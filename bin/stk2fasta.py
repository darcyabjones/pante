#!/usr/bin/env python3

import re
import sys
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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


def line_to_seq(line, regex):
    sline = regex.split(line.strip(), maxsplit=1)
    seq = SeqRecord(
        id=sline[0],
        name=sline[0],
        description=sline[0],
        seq=Seq(sline[1].replace("-", "").replace(".", "")),
    )

    return seq


def main():
    args = cli(sys.argv[0], sys.argv[1:])
    regex = re.compile(r"\s+")

    seen = set()
    seqs = []
    for line in args.infile:
        if line.startswith("#") or line.startswith("/") or len(line) == 1:
            continue

        seq = line_to_seq(line, regex)
        if seq.id not in seen:
            seqs.append(seq)

        seen.add(seq.id)

    SeqIO.write(seqs, args.outfile, "fasta")
    return


if __name__ == "__main__":
    main()
