#!/usr/bin/env python3

import sys
import argparse
from os.path import splitext

from Bio import SeqIO


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="""
        Splits a fasta file by sequence into N partitions of roughly
        the same size. Will not split one sequence into two.
        """
    )

    parser.add_argument(
        "infile",
        default=sys.stdin,
        type=argparse.FileType('r'),
        help="Input fasta files. Use '-' for stdin.",
    )

    parser.add_argument(
        "-p", "--prefix",
        default="out_",
        help="Output fasta file path. Default stdout.",
    )

    parser.add_argument(
        "-n", "--nchunks",
        default=8,
        type=int,
        help="The maximum number of partitions to separate into.",
    )

    return parser.parse_args(args)


def chunk_seqs(lengths, n=10, penalty=15):

    lengths = sorted(lengths.items(), key=lambda t: t[1], reverse=True)

    total_length = sum(l[1] for l in lengths)
    av_chunk_length = total_length / n

    chunks = [[] for i in range(n)]

    for i in range(len(lengths)):
        sid, length = lengths[i]
        dists = []
        for j in range(n):
            chunk_length = sum(t[1] for t in chunks[j]) + length
            diff = av_chunk_length - chunk_length
            # If projected length is greater than average, penalise it
            if diff < 0:
                diff *= -1
                diff *= penalty
            dists.append(diff)

        dist_min = min(dists)
        dist_min_idx = [
            j
            for j, d
            in enumerate(dists)
            if d == dist_min
        ][0]

        chunks[dist_min_idx].append(lengths[i])

    chunks.sort(key=lambda x: sum(z for y, z in x), reverse=True)
    return chunks


def main():
    args = cli(sys.argv[0], sys.argv[1:])
    seqs = SeqIO.to_dict(SeqIO.parse(args.infile, format="fasta"))
    lengths = {id_: len(seq.seq) for id_, seq in seqs.items()}

    ext = splitext(args.infile.name)[1]
    if ext == "" or ext is None:
        ext = ".fasta"

    chunks = chunk_seqs(
        lengths,
        n=args.nchunks,
        penalty=15
    )

    for i, chunk in enumerate(chunks):
        if len(chunk) == 0:
            continue

        these_seqs = [seqs[sid] for sid, length in chunk]
        outfile = f"{args.prefix}{i}{ext}"
        SeqIO.write(these_seqs, outfile, format="fasta")
    return


if __name__ == "__main__":
    main()
