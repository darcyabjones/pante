#!/usr/bin/env python3

import sys
import argparse
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
        "infile",
        type=argparse.FileType('r'),
        help="Input fasta file.",
    )

    parser.add_argument(
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help=("Output filtered fastas."),
    )

    parser.add_argument(
        "-c", "--coverage",
        default=0.95,
        type=float,
        help="The coverage of the smaller sequence required to exclude."
    )

    return parser.parse_args(args)


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


def coverage(left, right):
    intersection_begin = max([left.begin, right.begin])
    intersection_end = min([left.end, right.end])
    size = intersection_end - intersection_begin
    return size / min([left.length(), right.length()])


def main():
    args = cli(sys.argv[0], sys.argv[1:])
    regex = re.compile(r"^>?(?P<genome>[^:\s]+):(?P<seqid>[^:\s]+)"
                       r":(?P<start>\d+)-(?P<end>\d+)")

    itree = defaultdict(IntervalTree)
    seqs = SeqIO.to_dict(SeqIO.parse(args.infile, format="fasta"))

    for id_ in seqs.keys():
        genome, seqid, start, end = parse_id_as_interval(id_, regex)
        interval = Interval(start, end, data=id_)

        if interval in itree[(genome, seqid)]:
            continue

        envelops = itree[(genome, seqid)].envelop(interval)
        overlaps = itree[(genome, seqid)].overlap(interval)

        # If the interval completely overlaps one already in the tree
        # replace it. This would be covered by overlaps, by should be faster.
        if len(envelops) > 0:
            itree[(genome, seqid)].remove_overlap(start, end)
            itree[(genome, seqid)].add(interval)

        # If the interval partially overlaps one, or is completely contained
        # by an interval in the tree, we interrogate further.
        elif len(overlaps) > 0:
            to_remove = []
            add_to_tree = True

            for i in overlaps:
                cov_match = coverage(i, interval) > args.coverage

                # If the coverage of the shorter interval is above a threshold
                # and the interval already in the tree is the shorter one,
                # we flag it for replacement.
                if cov_match and i.length() < interval.length():
                    to_remove.append(i)

                # If the new interval was the shorter of the intervals.
                elif cov_match:
                    add_to_tree = False
                    break

            # We reached all the way through without discarding the new
            # interval.
            if add_to_tree:
                for i in to_remove:
                    itree[(genome, seqid)].remove(i)

                itree[(genome, seqid)].add(interval)

        # If it doesn't overlap the sequence at all.
        else:
            itree[(genome, seqid)].add(interval)

    for (genome, seqid), subitree in itree.items():
        SeqIO.write((keep[i.data] for i in subitree), args.outfile, format="fasta")

    return


if __name__ == "__main__":
    main()
