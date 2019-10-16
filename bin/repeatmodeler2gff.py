#!/usr/bin/env python3

import sys
import argparse

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
        "-s", "--source",
        default="RepeatModeler",
        type=str,
        help="The source to use.",
    )

    parser.add_argument(
        "-t", "--type",
        default="repeat_region",
        type=str,
        help="The type to use.",
    )

    parser.add_argument(
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help="Output fasta file path. Default stdout.",
    )

    return parser.parse_args(args)


def parse_block(handle, source, type_):

    note = []
    name = None
    species = []

    ids = []

    for line in handle:
        sline = line.strip()
        if line.startswith("//"):
            break

        if line.startswith("#=GF ID"):
            name = sline.split("ID", maxsplit=1)[-1].strip()

        elif line.startswith("#=GF DE"):
            note.append(sline.split("DE", maxsplit=1)[-1].strip())

        elif line.startswith("#=GF TP"):
            species.extend(
                sline.split("TP", maxsplit=1)[-1].strip().split(";")
            )

        elif not line.startswith("#"):
            ids.append(line.split(maxsplit=1)[0])

    out = []
    for id_ in ids:
        seqid, start, end, strand = id_to_loc(id_)

        if len(species) == 0:
            custom = {}
        else:
            custom = {"species": ":".join(species)}

        attributes = GFFAttributes(
            name=name,
            note=note,
            custom=custom,
        )

        record = GFFRecord(
            seqid=seqid,
            source=source,
            type=type_,
            start=start,
            end=end,
            score=None,
            strand=strand,
            attributes=attributes,
        )

        out.append(record)

    return out


def id_to_loc(string):
    seqid, loc = string.split(":", maxsplit=1)
    start_tmp, end_tmp = loc.split("-", maxsplit=1)
    start = int(start_tmp)
    end = int(end_tmp)

    if start <= end:
        return seqid, start - 1, end, Strand.PLUS
    else:
        return seqid, end - 1, start, Strand.MINUS


def main():
    args = cli(sys.argv[0], sys.argv[1:])

    block = parse_block(args.infile, args.source, args.type)
    while len(block) > 0:
        for record in block:
            print(record, file=args.outfile)

        block = parse_block(args.infile, args.source, args.type)
    return


if __name__ == "__main__":
    main()
