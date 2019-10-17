#!/usr/bin/env python3

import sys
import argparse
from collections import defaultdict

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

from gffpal.gff import GFF, Strand


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description=""" """
    )

    parser.add_argument(
        "ingff",
        type=argparse.FileType('r'),
        help=".",
    )

    parser.add_argument(
        "infasta",
        type=argparse.FileType('r'),
        help=".",
    )

    parser.add_argument(
        "-s", "--source",
        type=str,
        default="OcculterCut",
        help="The source.",
    )

    parser.add_argument(
        "-t", "--type",
        type=str,
        default="region",
        help="The type.",
    )

    parser.add_argument(
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help="Output GFF3 file. Default stdout.",
    )

    return parser.parse_args(args)


def gff_to_seqfeature(gff):
    loc = FeatureLocation(gff.start, gff.end, 1)
    return SeqFeature(loc)


def reverse_complement(string):
    out = []
    for char in reversed(string):
        if char == "A":
            out.append("T")
        elif char == "T":
            out.append("A")
        elif char == "G":
            out.append("C")
        elif char == "C":
            out.append("G")
        else:
            out.append("N")

    return "".join(out)


def count_frequencies(seq):

    mono = defaultdict(int)
    di = defaultdict(int)

    last = None
    for base in str(seq.seq).upper():
        mono[base] += 1
        if last is not None:
            dn = last + base
            rc_dn = reverse_complement(dn)
            di[min([dn, rc_dn])] += 1

        last = base

    out = {}
    mono_total = sum(mono.values())
    out["at"] = (mono["A"] + mono["T"]) / mono_total
    out["gc"] = (mono["G"] + mono["C"]) / mono_total

    di_total = sum(di.values())
    for dn, count in di.items():
        if "N" in dn:
            continue

        name = (dn[0] + "p" + dn[1]).lower()
        out[name] = count / di_total

    return out


def main():
    args = cli(sys.argv[0], sys.argv[1:])

    gff = GFF.parse(args.ingff)
    seqs = SeqIO.to_dict(SeqIO.parse(args.infasta, format="fasta"))

    counter = 1
    for region in gff:
        region.type = args.type
        region.source = args.source
        region.strand = Strand.UNSTRANDED

        region.attributes.id = f"{args.type}{counter}"

        sf = gff_to_seqfeature(region)
        seq = sf.extract(seqs[region.seqid])

        base_counts = count_frequencies(seq)
        region.attributes.custom = base_counts

        counter += 1

        print(region, file=args.outfile)
    return


if __name__ == "__main__":
    main()
