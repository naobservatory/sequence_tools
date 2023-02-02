#!/usr/bin/env python3

import re
import sys
import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser

def start():
    parser = argparse.ArgumentParser(
        description='Take fasta marked up with alignment information and '
        'display them horizontally aligned')
    parser.add_argument('--input',
                        metavar='foo.aligned.fasta',
                        default="/dev/stdin")
    parser.add_argument(
        '--sort',
        action=argparse.BooleanOptionalAction,
        default=True,
        help='Sort lines before printing.  Use when order in the input file '
        'is not meaningful.')
    parser.add_argument(
        "contig",
        help='The contig the input was aligned to, for displaying as a header')
    run(parser.parse_args())

def run(args):
    reads = []
    with open(args.input) as inf:
        for seqid, seq in SimpleFastaParser(inf):
            pos, = re.findall("pos=([-0-9]+)", seqid)
            reads.append((int(pos), seq))

    if not reads:
        raise Exception("no reads: supply on stdin or set --input")

    if args.sort:
        reads.sort()

    contig_offset = -sorted(reads)[0][0]
    if contig_offset < 0:
        contig_offset = 0

    print("%s%s" % (" "*contig_offset, args.contig))

    for contig_pos, read in reads:
        print("%s%s" % (
            " "*(contig_pos + contig_offset), read))

if __name__ == "__main__":
    start()
