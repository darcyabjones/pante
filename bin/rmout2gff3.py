#!/usr/bin/env python3

import re
import sys

import argparse

from gffpal.gff import GFFRecord, Strand
from gffpal.attributes import Target, GFFAttributes


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description=""" Converts a space sepatated out file from repeatmasker
        to a GFF3.
        """
    )

    parser.add_argument(
        "infile",
        type=argparse.FileType('r'),
        help="Input '.out' file.",
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
        default="RepeatMasker",
        help="A tab-separated table mapping query IDs to gff3 attributes.",
    )

    parser.add_argument(
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help="Output GFF3 file. Default stdout.",
    )

    return parser.parse_args(args)


class RMOut(object):

    def __init__(
        self,
        score: int,
        perc_divergence: float,
        perc_deletions: float,
        perc_insertions: float,
        query: str,
        qstart: int,
        qend: int,
        qremaining: int,
        strand: str,
        family: str,
        kind: str,
        tstart: int,
        tend: int,
        tremaining: int,
        index: int,
        better_hit: bool,
    ):

        self.score = score
        self.perc_divergence = perc_divergence
        self.perc_deletions = perc_divergence
        self.perc_insertions = perc_insertions
        self.query = query
        self.qstart = qstart
        self.qend = qend
        self.qremaining = qremaining
        self.strand = strand
        self.family = family
        self.kind = kind
        self.tstart = tstart
        self.tend = tend
        self.tremaining = tremaining
        self.index = index
        self.better_hit = better_hit
        return

    @classmethod
    def parse(cls, line):
        sline = re.split(r"\s+", line.strip())

        sw = int(sline[0])
        perc_divergence = float(sline[1])
        perc_deletions = float(sline[2])
        perc_insertions = float(sline[3])
        query = str(sline[4])
        qstart = int(sline[5]) - 1
        qend = int(sline[6])
        qremaining = int(sline[7].strip("()"))
        strand = "-" if sline[8] == "C" else "+"
        family = sline[9]
        kind = sline[10]

        if strand == "+":
            tstart = int(sline[11]) - 1
            tend = int(sline[12])
            tremaining = int(sline[13].strip("()"))
        else:
            tstart = int(sline[13]) - 1
            tend = int(sline[12])
            tremaining = int(sline[11].strip("()"))

        index = int(sline[14])

        if len(sline) > 15:
            better_hit = sline[15] == "*"
        else:
            better_hit = False

        return cls(
            sw,
            perc_divergence,
            perc_deletions,
            perc_insertions,
            query,
            qstart,
            qend,
            qremaining,
            strand,
            family,
            kind,
            tstart,
            tend,
            tremaining,
            index,
            better_hit,
        )

    @classmethod
    def from_file(cls, handle):
        started = False

        for line in handle:
            sline = line.strip()
            if not started:
                started = sline == ""
                continue

            yield cls.parse(sline)
        return

    def as_gff3(self, source="RepeatMasker"):

        custom = {
            "smith_waterman_score": str(self.score),
            "percent_divergence": str(self.perc_divergence),
            "percent_deletions": str(self.perc_deletions),
            "percent_insertions": str(self.perc_insertions),
            "family_consensus_length": str(self.tend + self.tremaining),
        }

        ontology_terms = []

        if self.better_hit:
            custom["has_better_overlapping_hit"] = "true"

        if self.kind == "Simple_repeat":
            repeat_unit = re.match(
                r"\((?P<rep>.+)\)n",
                self.family
            ).group("rep")

            custom["repeat_unit"] = repeat_unit

            if len(repeat_unit) == 1:
                type_ = "monomeric_repeat"
                ontology_terms.append("SO:0001934")
            elif len(repeat_unit) < 10:
                type_ = "microsatellite"
                ontology_terms.append("SO:0000289")
            else:
                type_ = "minisatellite"
                ontology_terms.append("SO:0000643")

        elif self.kind == "Low_complexity":
            type_ = "low_complexity_region"
            ontology_terms.append("SO:0001005")

        elif self.kind == "Other/DNA_virus":
            type_ = "viral_sequence"
            ontology_terms.append("SO:0001041")
            ontology_terms.append("SO:0000347")
            custom["repeat_family"] = self.kind

        elif self.kind.startswith("snRNA"):
            type_ = "snRNA_gene"
            ontology_terms.append("SO:0001268")
            ontology_terms.append("SO:0000347")
            custom["repeat_family"] = self.kind

        elif self.kind.startswith("tRNA"):
            type_ = "tRNA_gene"
            ontology_terms.append("SO:0001272")
            ontology_terms.append("SO:0000347")
            custom["repeat_family"] = self.kind

        elif self.kind.startswith("rRNA"):
            type_ = "rRNA_gene"
            ontology_terms.append("SO:0000347")
            ontology_terms.append("SO:0001637")
            custom["repeat_family"] = self.kind

        elif self.kind.startswith("scRNA"):
            type_ = "scRNA_gene"
            ontology_terms.append("SO:0000347")
            ontology_terms.append("SO:0001266")
            custom["repeat_family"] = self.kind

        elif self.kind.startswith("Segmental"):
            type_ = "duplication"
            ontology_terms.append("SO:0000347")
            ontology_terms.append("SO:1000035")
            custom["repeat_family"] = self.kind

        elif self.kind.startswith("Satellite"):
            type_ = "satellite_DNA"
            ontology_terms.append("SO:0000347")
            ontology_terms.append("SO:0000005")
            custom["repeat_family"] = self.kind

        elif self.kind.startswith("Retrotransposon"):
            type_ = "retrotransposon"
            ontology_terms.append("SO:0000347")
            ontology_terms.append("SO:0000180")
            custom["repeat_family"] = self.kind

        elif self.kind.startswith("DNA"):
            type_ = "DNA_transposon"
            ontology_terms.append("SO:0000182")
            ontology_terms.append("SO:0000347")
            custom["repeat_family"] = self.kind

        elif self.kind.startswith("LTR"):
            type_ = "LTR_retrotransposon"
            ontology_terms.append("SO:0000186")
            ontology_terms.append("SO:0000347")
            custom["repeat_family"] = self.kind

        elif self.kind.startswith("LINE"):
            type_ = "LINE_element"
            ontology_terms.append("SO:0000194")
            ontology_terms.append("SO:0000347")
            custom["repeat_family"] = self.kind

        elif self.kind.startswith("SINE"):
            type_ = "SINE_element"
            ontology_terms.append("SO:0000206")
            ontology_terms.append("SO:0000347")
            custom["repeat_family"] = self.kind

        else:
            type_ = "nucleotide_match"
            ontology_terms.append("SO:0000347")
            custom["repeat_family"] = self.kind

        if "helitron" in self.kind.lower():
            type_ = "helitron"
            ontology_terms.append("SO:0000544")

        if ("maverick" in self.kind.lower()
                or "polinton" in self.kind.lower()):
            type_ = "polinton"
            ontology_terms.append("SO:0001170")

        if "mite" in self.kind.lower():
            type_ = "MITE"
            ontology_terms.append("SO:0000338")

        if (self.kind.lower().endswith("/P")
                or self.kind.lower().endswith("P-Fungi")):
            type_ = "p_element"
            ontology_terms.append("SO:0001535")

        attributes = GFFAttributes(
            name=self.family,
            target=Target(self.family, self.tstart, self.tend),
            ontology_term=ontology_terms,
            custom=custom,
        )

        if type_ in {
                "monomeric_repeat", "microsatellite", "minisatellite",
                "low_complexity_region", "satellite_DNA", "duplication"
        }:
            strand = Strand.UNSTRANDED
        else:
            strand = Strand.parse(self.strand)

        record = GFFRecord(
            seqid=self.query,
            source=source,
            type=type_,
            start=self.qstart,
            end=self.qend,
            score=self.perc_divergence,
            strand=strand,
            attributes=attributes,
        )
        return record


def main():
    args = cli(sys.argv[0], sys.argv[1:])
    lines = RMOut.from_file(args.infile)

    for line in lines:
        if args.best and line.better_hit:
            continue

        print(line.as_gff3(source=args.source), file=args.outfile)
    return


if __name__ == "__main__":
    main()
