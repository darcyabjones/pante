#!/usr/bin/env python3

import sys
import argparse
import os
import re

from collections import defaultdict

from Bio import SeqIO

from intervaltree import Interval, IntervalTree


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description=""" Extracts all sequences from a fasta msa file into a
        stockholm alignment file.
        """
    )

    parser.add_argument(
        "table",
        type=argparse.FileType('r'),
        help="Input table file.",
    )

    parser.add_argument(
        "infiles",
        type=argparse.FileType('r'),
        nargs="+",
        help="Input fasta files.",
    )

    parser.add_argument(
        "-s", "--outsel",
        default="selected_families",
        type=str,
        help=(
            "Output directory for families that pass the tests. "
            "Default 'selected_families'."
        ),
    )

    parser.add_argument(
        "-e", "--outexcl",
        default="excluded_families",
        type=str,
        help=(
            "Output directory for families that don't pass the tests. "
            "Default 'excluded_families'."
        ),
    )

    parser.add_argument(
        "-t", "--outtab",
        default="filtered_families.tsv",
        type=argparse.FileType('w'),
        help=("Output table summarising results. "
              "Default 'filtered_families.tsv'."),
    )

    parser.add_argument(
        "-m", "--min-intra",
        default=5,
        type=int,
        help="The minimum number of repeat family loci within a single genome."
    )

    parser.add_argument(
        "-c", "--min-inter",
        default=0.05,
        type=float,
        help="The minimum proportion of isolates that should have this repeat."
    )

    parser.add_argument(
        "-n", "--size",
        default=None,
        type=int,
        help=("The number of isolates in the population. "
              "By default will estimate from table.")
    )

    return parser.parse_args(args)


class Table(object):

    def __init__(
        self,
        rtype,
        cluster,
        length,
        psim,
        strand,
        cigar,
        member,
        centroid
    ):
        self.rtype = rtype
        self.cluster = cluster
        self.length = length
        self.psim = psim
        self.strand = strand
        self.cigar = cigar
        self.member = member
        self.centroid = centroid
        return

    @classmethod
    def parse(cls, handle):
        columns = ["rtype", "cluster", "length", "psim", "strand",
                   "mis1", "mis2", "cigar", "member", "centroid"]

        for line in handle:
            sline = line.strip().split("\t")
            assert len(sline) == len(columns)

            row = dict(zip(columns, sline))
            yield cls(
                rtype=row["rtype"],
                cluster=row["cluster"],
                length=int(row["length"]),
                psim=None if row["psim"] == "*" else float(row["psim"]),
                strand=row["strand"],
                cigar=row["cigar"],
                member=row["member"],
                centroid=None if row["centroid"] == "*" else row["centroid"],
            )

        return


def parse_id_as_interval(id_string, regex):
    """ The fasta ids contain the locus information. """

    match = regex.match(id_string)
    genome = match.group("genome")
    seqid = match.group("seqid")
    start_tmp = int(match.group("start"))
    end_tmp = int(match.group("end"))

    start = min([start_tmp, end_tmp])
    end = max([start_tmp, end_tmp])
    del start_tmp
    del end_tmp

    return (genome, seqid, start, end)


def get_unique_loci(intervals):

    grouped_intervals = defaultdict(list)
    for genome, seqid, start, end in intervals:
        grouped_intervals[(genome, seqid)].append(Interval(start, end))

    unique_loci = list()
    for (genome, seqid), intvls in grouped_intervals.items():
        itree = IntervalTree(intvls)
        itree.merge_overlaps()
        for intvl in itree:
            unique_loci.append((genome, seqid, intvl.begin, intvl.end))

    return unique_loci


def find_intra_counts(loci):

    count = defaultdict(int)
    for genome, seqid, start, end in loci:
        count[genome] += 1

    return count


def filter_intra_counts(intra_counts, min_count=2):
    return {g: c for g, c in intra_counts.items() if c >= min_count}


def find_inter_count(intra_counts):

    count = 0
    for genome, intra_count in intra_counts.items():
        count += 1

    return count


def get_count_statistics(intervals, min_intra=2):

    unique_loci = get_unique_loci(intervals)
    intra_counts = find_intra_counts(unique_loci)

    filtered_intra_counts = filter_intra_counts(
        intra_counts,
        min_intra
    )

    inter_count = find_inter_count(filtered_intra_counts)
    return unique_loci, intra_counts, inter_count


def get_size_from_table(table, regex):
    genomes = set()

    for row in table:
        g, c, s, e = parse_id_as_interval(row.member, regex)
        genomes.add(g)

    return len(genomes)


def write_outtab_header(outfile):
    line = [
        "num",
        "centroid",
        "nunique",
        "inter",
        "infreq",
        "genome",
        "intra",
    ]

    print("\t".join(line), file=outfile)
    return


def write_outtab(
    outfile,
    num,
    centroid,
    nunique,
    inter,
    infreq,
    genome,
    intra
):
    line = [num, centroid, str(nunique), str(inter),
            str(infreq), genome, str(intra)]

    print("\t".join(line), file=outfile)
    return


def get_id_to_table_member(table):
    out = dict()

    for line in table:
        if line.rtype != "H":
            continue

        out[line.member] = line

    return out


def get_id_to_centroid(table):
    out = dict()

    for line in table:
        if line.centroid is None:
            centroid = line.member
        else:
            centroid = line.centroid

        out[line.member] = centroid
    return out


def get_id_to_num(table):
    out = dict()

    for line in table:
        print(line.rtype, line.member, line.centroid, line.cluster)
        if line.member is not None:
            out[line.member] = line.cluster
        elif line.centroid is not None:
            out[line.centroid] = line.cluster

    return out


def main():
    args = cli(sys.argv[0], sys.argv[1:])
    regex = re.compile(r"^>?(?P<genome>[^:\s]+):(?P<seqid>[^:\s]+)"
                       r":(?P<start>\d+)-(?P<end>\d+)")

    table = list(Table.parse(args.table))
    id_to_table_member = get_id_to_table_member(table)
    id_to_num = get_id_to_num(table)
    id_to_centroid = get_id_to_centroid(table)

    if args.size is None:
        size = get_size_from_table(table, regex)
    else:
        size = args.size

    write_outtab_header(args.outtab)

    os.mkdir(args.outsel)
    os.mkdir(args.outexcl)

    for infile in args.infiles:

        seqs = SeqIO.to_dict(SeqIO.parse(
            infile,
            format="fasta"
        ))

        cluster_nums = set()
        for id_ in seqs.keys():
            num = id_to_num.get(id_, None)

            if num is None:
                raise ValueError(f"Missing in id_to_num {id_}")

            cluster_nums.add(num)

        cluster_num = list(cluster_nums)[0]

        cluster_centroids = set()
        for id_ in seqs.keys():
            centroid = id_to_centroid.get(id_, None)

            if centroid is None:
                raise ValueError(f"Missing in id_to_centroid {id_}")
            cluster_centroids.add(centroid)

        cluster_centroid = list(cluster_centroids)[0]

        intervals = {
            id: parse_id_as_interval(id, regex)
            for id
            in seqs.keys()
        }

        unique_loci, intra_counts, inter_count = get_count_statistics(
            intervals.values(),
            min_intra=args.min_intra
        )

        inter_freq = inter_count / size

        for genome, intra_count in intra_counts.items():
            write_outtab(
                args.outtab,
                cluster_num,
                cluster_centroid,
                len(unique_loci),
                inter_count,
                inter_freq,
                genome,
                intra_count,
            )

        out_seqs = []
        for id_, seq in seqs.items():
            tab = id_to_table_member.get(id_, None)

            if tab is None:
                strand = "+"
            elif tab.strand == "-":
                strand = "-"
            else:
                strand = "+"

            if strand == "-":
                seq = seq.reverse_complement()
                # Flip the strand.
                g, c, s, e = intervals[id_]
                new_id = f"{g}:{c}:{e}-{s}"

                seq.id = new_id
                seq.name = new_id
                seq.description = new_id

            out_seqs.append(seq)

        if inter_freq >= args.min_inter:
            outdir = args.outsel
        else:
            outdir = args.outexcl

        outpath = os.path.join(outdir, os.path.split(infile.name)[-1])
        SeqIO.write(out_seqs, outpath, format="fasta")

    return


if __name__ == "__main__":
    main()
