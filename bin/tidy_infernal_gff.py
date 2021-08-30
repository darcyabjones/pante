#!/usr/bin/env python3

import sys
import argparse
from collections import defaultdict

from gffpal.gff import GFFRecord
from gffpal.attributes import Target, GFFAttributes


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description=""" """
    )

    parser.add_argument(
        "infile",
        type=argparse.FileType('r'),
        help=".",
    )

    parser.add_argument(
        "-g", "--go",
        type=argparse.FileType('r'),
        default=None,
        help="rfam2go",
    )

    parser.add_argument(
        "-b", "--best",
        action="store_true",
        default=False,
        help="Exclude matches that overlap a better scoring match.",
    )

    parser.add_argument(
        "-s", "--source",
        type=str,
        default="infernal",
        help="The source.",
    )

    parser.add_argument(
        "-t", "--type",
        type=str,
        default="ncRNA_gene",
        help="The type to use for matches.",
    )

    parser.add_argument(
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help="Output GFF3 file. Default stdout.",
    )

    return parser.parse_args(args)


def parse_rfam2go(handle):
    out = defaultdict(list)
    for line in handle:
        line = line.strip()
        if line.startswith("!"):
            continue
        
        rfam, go = line.strip().split(">", maxsplit=1)
        rfam_acc, rfam_name = rfam.strip().split()
        go_acc, go_name = go.strip().split(";", maxsplit=1)

        rfam_acc_stripped = rfam_acc.strip().split(":", maxsplit=1)[-1].strip()
        out[rfam_acc_stripped].extend([go_acc.strip(), go_name.strip()])

    return out


def main():

    args = cli(sys.argv[0], sys.argv[1:])

    if args.go is not None:
        go = parse_rfam2go(args.go)
    else:
        go = {}

    for line in args.infile:
        if line.startswith("#"):
            continue

        record = GFFRecord.parse(line)
        attrs = record.attributes

        if args.best and attrs.custom["olp"] == "=":
            continue

        name = record.type
        dbxrefs = ["Rfam:" + attrs.custom["mdlaccn"], "Rfam:" + name]
        if "clan" in attrs.custom:
            dbxrefs.append("RfamClan:" + attrs.custom["clan"])

        ontology_terms = go.get(attrs.custom["mdlaccn"], [])
        notes = [attrs.custom["desc"]]

        target = Target(
            attrs.custom["mdlaccn"],
            int(attrs.custom["mdlfrom"]),
            int(attrs.custom["mdlto"]),
        )

        custom = {
            "evalue": attrs.custom["evalue"],
            "model_type": attrs.custom["mdl"],
            "gc": attrs.custom["gc"],
            "bias": attrs.custom["bias"],
            "bitscore": record.score,
        }

        if attrs.custom["trunc"] == "yes":
            custom["truncated_match"] = "true"

        if attrs.custom["olp"] == "=":
            custom["overlap_with_better_score"] = "true"

        record.source = args.source
        record.type = args.type
        record.score = float(custom["evalue"])
        record.attributes = GFFAttributes(
            name=name,
            dbxref=dbxrefs,
            target=target,
            note=notes,
            ontology_term=ontology_terms,
            custom=custom,
        )

        print(record, file=args.outfile)


if __name__ == "__main__":
    main()
