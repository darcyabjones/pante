#!/usr/bin/env python3

import sys
import argparse

from Bio import AlignIO

from Bio.Alphabet import Gapped, SingleLetterAlphabet

from gffpal.gff import GFFRecord, Strand
from gffpal.attributes import GFFAttributes


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
    alignments = AlignIO.parse(
        args.infile,
        format="stockholm",
        alphabet=Gapped(SingleLetterAlphabet(), "-")
    )

    i = 1
    for alignment in alignments:
        name = f"fam{i}"
        for sequence in alignment:
            seqid, loc = sequence.id.split(":", maxsplit=1)
            start_tmp, end_tmp = loc.split("-", maxsplit=1)
            start = int(start_tmp)
            end = int(end_tmp)

            if start <= end:
                strand = Strand.PLUS
            else:
                strand = Strand.MINUS
                tmp = end
                end = start
                start = tmp
                del tmp

            attributes = GFFAttributes(name=name)
            record = GFFRecord(
                seqid=seqid,
                source="RepeatModeler",
                type="repeat_region",
                start=start,
                end=end,
                score=None,
                strand=strand,
                attributes=attributes,
            )

            print(record, file=args.outfile)
        i += 1
    return


if __name__ == "__main__":
    main()
