#!/usr/bin/env python3

import sys
import argparse

from Bio import SeqIO

from gffpal.gff import GFFRecord, Strand
from gffpal.attributes import GFFAttributes


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description=""" Converts a tab-separated blast-like file to a GFF3.

        The table should have a header, and column names should match
        mmseqs labels.
        """
    )

    parser.add_argument(
        "genome",
        type=argparse.FileType('r'),
        help="The source.",
    )

    parser.add_argument(
        "infile",
        type=argparse.FileType('r'),
        help="Input fasta file.",
    )

    parser.add_argument(
        "-s", "--source",
        type=str,
        default="mitefinder",
        help="The source.",
    )

    parser.add_argument(
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help="Output GFF3 file. Default stdout.",
    )

    return parser.parse_args(args)


def split_desc(seq, seqids):
    split_id = seq.id.split("|")
    seqid = seqids[int(split_id[1]) - 1]
    lborder_start = int(split_id[2])
    lborder_end = int(split_id[3])
    rborder_start = int(split_id[4])
    rborder_end = int(split_id[5])
    score = float(split_id[-1].split(":", maxsplit=1)[1])
    return seqid, lborder_start, lborder_end, rborder_start, rborder_end, score


def get_region_feature(i, seqid, left, right, score):

    start = min(left + right)
    end = max(left + right)

    region_id = f"repeat_region{i}"
    attributes = GFFAttributes(
        id=region_id,
        ontology_term=["SO:0000657", "SO:repeat_region"],
        custom={"mitefinder_score": score},
    )

    record = GFFRecord(
        seqid=seqid,
        source="MiteFinderII",
        type="repeat_region",
        start=start,
        end=end,
        score=score,
        strand=Strand.UNKNOWN,
        attributes=attributes,
    )

    return record


def get_tir_feature(i, seqid, pos):
    start = min(pos)
    end = max(pos)
    region_id = f"repeat_region{i}"
    attributes = GFFAttributes(
        parent=[region_id],
        ontology_term=["SO:0000481", "SO:terminal_inverted_repeat"]
    )

    record = GFFRecord(
        seqid=seqid,
        source="MiteFinderII",
        type="terminal_inverted_repeat",
        start=start,
        end=end,
        score=None,
        strand=Strand.UNKNOWN,
        attributes=attributes,
    )

    return record


def get_mite_feature(i, seqid, left, right):
    start = max(left)
    end = min(right)
    assert start < end

    region_id = f"repeat_region{i}"

    attributes = GFFAttributes(
        parent=[region_id],
        ontology_term=["SO:0000338", "SO:MITE"]
    )

    record = GFFRecord(
        seqid=seqid,
        source="MiteFinderII",
        type="MITE",
        start=start,
        end=end,
        score=None,
        strand=Strand.UNKNOWN,
        attributes=attributes,
    )
    return record


def main():
    args = cli(sys.argv[0], sys.argv[1:])
    seqs = SeqIO.parse(args.infile, format="fasta")
    genome = SeqIO.parse(args.genome, format="fasta")
    seqids = [s.id for s in genome]

    i = 1
    for seq in seqs:
        (seqid, lborder_start, lborder_end,
         rborder_start, rborder_end, score) = split_desc(seq, seqids)

        region = get_region_feature(
            i,
            seqid,
            [lborder_start, lborder_end],
            [rborder_start, rborder_end],
            score
        )

        ltir = get_tir_feature(i, seqid, [lborder_start, lborder_end])
        rtir = get_tir_feature(i, seqid, [rborder_start, rborder_end])

        mid = get_mite_feature(
            i,
            seqid,
            [lborder_start, lborder_end],
            [rborder_start, rborder_end],
        )

        print(region, file=args.outfile)
        print(ltir, file=args.outfile)
        print(mid, file=args.outfile)
        print(rtir, file=args.outfile)
        i += 1
    return


if __name__ == "__main__":
    main()
