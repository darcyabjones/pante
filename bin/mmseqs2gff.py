#!/usr/bin/env python3

import sys
import csv
import re
import argparse

from gffpal.gff import GFFRecord, Strand
from gffpal.attributes import Target, Gap, GapCode, GapElement, GFFAttributes


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description=""" Converts a tab-separated blast-like file to a GFF3.

        The table should have a header, and column names should match
        mmseqs labels.
        """
    )

    parser.add_argument(
        "infile",
        type=str,
        help="Input fasta files.",
    )

    parser.add_argument(
        "-a", "--attributes",
        type=str,
        default=None,
        help="A tab-separated table mapping query IDs to gff3 attributes.",
    )

    parser.add_argument(
        "-s", "--source",
        type=str,
        default="pante",
        help="A tab-separated table mapping query IDs to gff3 attributes.",
    )

    parser.add_argument(
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help="Output GFF3 file. Default stdout.",
    )

    return parser.parse_args(args)


def parse_cigar(cigar):
    output = []

    parts = re.findall(r"\d+[MIDFR]", cigar)
    for part in parts:
        num = int(part[:-1])
        code = GapCode.parse(part[-1])
        output.append(GapElement(code, num))
    return Gap(output)


def parse_external_attrs(handle):
    ext_attrs = {}
    reader = csv.DictReader(handle, delimiter="\t")
    for line in reader:
        if line["Name"] == ".":
            line["Name"] = None

        if line["Alias"] == ".":
            line["Alias"] = []
        else:
            line["Alias"] = line["Alias"].split(",")

        if line["Note"] == ".":
            line["Note"] = []
        else:
            line["Note"] = [line["Note"]]

        if line["Ontology_term"] == ".":
            line["Ontology_term"] = []
        else:
            line["Ontology_term"] = line["Ontology_term"].split(",")

        if line["Dbxref"] == ".":
            line["Dbxref"] = []
        else:
            line["Dbxref"] = line["Dbxref"].split(",")

        ext_attrs[line["ID"]] = line
    return ext_attrs


def dict_to_record(row, source="gffpal", ext_attributes=None):

    if "#target" in row:
        seqid = row["#target"]
    elif "target" in row:
        seqid = row["target"]
    else:
        raise ValueError("Input file doesn't have a 'target' column.")

    start = int(row["tstart"])
    end = int(row["tend"])

    if start <= end:
        strand = Strand.PLUS
    else:
        tmp = start
        start = end
        end = tmp
        del tmp
        strand = Strand.MINUS

    custom_attrs = {
        "pident": row["pident"],
        "alnlen": row["alnlen"],
        "score": row["raw"],
        "bitscore": row["bits"],
        "gapopen": row["gapopen"],
        "query_coverage": row["qcov"],
        "evalue": row["evalue"],
    }

    target = Target(row["query"], int(row["qstart"]), int(row["qend"]))
    gap = parse_cigar(row["cigar"])

    if ext_attributes is None:
        attributes = GFFAttributes(target=target, custom=custom_attrs, gap=gap)
    else:
        attributes = GFFAttributes(
            target=target,
            custom=custom_attrs,
            gap=gap,
            name=ext_attributes[row["query"]]["Name"],
            alias=ext_attributes[row["query"]]["Alias"],
            dbxref=ext_attributes[row["query"]]["Dbxref"],
            ontology_term=ext_attributes[row["query"]]["Ontology_term"],
            note=ext_attributes[row["query"]]["Note"],
        )

    record = GFFRecord(
        seqid=seqid,
        source=source,
        type="nucleotide_to_protein_match",
        start=start,
        end=end,
        score=float(row["evalue"]),
        strand=strand,
        attributes=attributes,
    )
    return record


def main():
    args = cli(sys.argv[0], sys.argv[1:])

    try:
        infile = open(args.infile, "r", newline="")

        reader = csv.DictReader(infile, delimiter="\t")

        if args.attributes is None:
            ext_attrs = None
        else:
            attributes = open(args.attributes, "r", newline="")
            ext_attrs = parse_external_attrs(attributes)

        for line in reader:
            print(
                dict_to_record(
                    line,
                    source=args.source,
                    ext_attributes=ext_attrs,
                ),
                file=args.outfile
            )
    finally:
        infile.close()

        if args.attributes is not None:
            attributes.close()

    return


if __name__ == "__main__":
    main()
