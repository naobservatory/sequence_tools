#!/usr/bin/env python3

import os
import re
import sys
import numpy as np
import argparse
import itertools

from Bio import Align
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord

from util import die, get_columns, interpret_sequence_argument

COLOR_RED = '\x1b[1;31m'
COLOR_GREEN = '\x1b[1;32m'
COLOR_YELLOW = '\x1b[1;33m'
COLOR_END = '\x1b[0m'

def wrap(s, columns):
    """Dead-simple wrapper."""
    line = []
    for c in s:
        line.append(c)
        if len(line) == columns:
            yield ''.join(line)
            line = []
    if line:
        yield ''.join(line)

def collapse_subs(alignment, max_dist):
    alignment_str = str(alignment)
    if "target" in alignment_str:
        # this needs relatively recent biopython
        rows = str(alignment).split('\n')
        targets = [row for row in rows if row.startswith("target")]
        queries = [row for row in rows if row.startswith("query")]

        def extract_seq(row):
            pieces = row.split()
            while pieces[-1].isdigit():
                pieces.pop()
            return pieces[-1]

        seq1 = "".join(extract_seq(row) for row in targets)
        seq2 = "".join(extract_seq(row) for row in queries)
    else:
        seq1, _, seq2, _ = alignment_str.split('\n')

    out1 = []
    out2 = []

    if len(seq1) != len(seq2):
        raise Exception(
            'sequences different length after alignment: %s vs %s' % (
                len(seq1), len(seq2)))

    for c1, c2 in zip(seq1, seq2):
        for i in range(1, min(max_dist, len(out1), len(out2))):
            if c1 == '-' and out2 and out2[-i] == '-':
                out2.pop(-i)
                out2.append(c2)
                break
            elif c2 == '-' and out1 and out1[-i] == '-':
                out1.pop(-i)
                out1.append(c1)
                break
        else:
            out1.append(c1)
            out2.append(c2)

    return ''.join(out1), ''.join(out2)

def color_mismatches(
        seq1_line, seq2_line,
        ins_color=COLOR_GREEN,
        del_color=COLOR_RED,
        sub_color=COLOR_YELLOW,
        match_color=COLOR_END):

    out1 = []
    out2 = []

    ins_color = COLOR_END + ins_color
    del_color = COLOR_END + del_color
    sub_color = COLOR_END + sub_color
    if match_color != COLOR_END:
        match_color = COLOR_END + match_color
        out1.append(match_color)
        out2.append(match_color)

    for c1, c2 in zip(seq1_line, seq2_line):
        c1_color = None
        c2_color = None
        if c1 == '-':
            c2_color = ins_color
        elif c2 == '-':
            c1_color = del_color
        elif c1 != c2:
            c1_color = sub_color
            c2_color = sub_color

        if c1_color:
            out1.append(c1_color)
        out1.append(c1)
        if c1_color:
            out1.append(match_color)

        if c2_color:
            out2.append(c2_color)
        out2.append(c2)
        if c2_color:
            out2.append(match_color)
    if match_color != COLOR_END:
        out1.append(COLOR_END)
        out2.append(COLOR_END)

    return ''.join(out1), ''.join(out2)

color_regex  =re.compile(r'\x1b\[[01](;3[123])?m')
def rmcolor(line):
    return color_regex.sub('', line)

def run(args):
    recs1 = interpret_sequence_argument(args.in1)
    recs2 = interpret_sequence_argument(args.in2)

    for rec1, rec2 in itertools.product(recs1, recs2):
        align_and_print(rec1, rec2, args)

def count_bases_in_printable_line(line):
    return len(rmcolor(line).replace("-", ""))

def print_moving_average_tsv(seq1_aligned, seq2_aligned, args):

    headers = ["position", "match", "average"]
    rows = []
    for i, (b1, b2) in enumerate(zip(seq1_aligned, seq2_aligned)):
        rows.append([i, 1 if b1 == b2 else 0, None])

    # set the weights to a gaussian centered at zero
    weights = {}
    normal_distribution_mean = 0
    standard_deviation = args.moving_average_width
    for i in range(-100, 101):
        weights[i] = (np.pi * standard_deviation) * np.exp(
            -0.5*((i-normal_distribution_mean)/standard_deviation)**2)

    for i in range(len(rows)):
        t = 0
        s = 0
        for j, w in weights.items():
            if 0 <= i + j < len(rows):
                t += w
                s += rows[i+j][1]*w

        rows[i][2] = s/t

    print("\t".join(headers))
    for position, match, average in rows:
        print("%s\t%s\t%.4f" % (position, match, average))

def trim_seqs(seq_a, seq_b, start_pos_a):
    progress_a = 0
    trimmed_a = []
    trimmed_b = []
    offset_a = 0
    offset_b = 0
    for a, b in zip(seq_a, seq_b):
        if a != "-":
            progress_a += 1
        if progress_a > start_pos_a:
            trimmed_a.append(a)
            trimmed_b.append(b)
        else:
            if a != "-":
                offset_a += 1
            if b != "-":
                offset_b += 1

    return "".join(trimmed_a), "".join(trimmed_b), offset_a, offset_b

def align_and_print(rec1, rec2, args):
    seq1, seq2 = rec1.seq, rec2.seq

    aligner = Align.PairwiseAligner()
    # These are the scoring settings porechop uses by default.
    # https://github.com/rrwick/Porechop/blob/master/porechop/porechop.py#L145
    aligner.end_gap_score = 0
    aligner.match_score = 3
    aligner.mismatch_score = -6
    aligner.internal_open_gap_score = -5
    aligner.internal_extend_gap_score = -2

    try:
        alignment = aligner.align(seq1, seq2)[0]
    except Exception:
        print(rec1.description)
        print(rec2.description)
        raise

    if alignment.score / min(len(seq1), len(seq2)) < (args.min_score/100):
        # print('no match: score=%s [%s, %s]' % (
        #     alignment.score, len(seq1), len(seq2)))
        return

    seq1_aligned, seq2_aligned = collapse_subs(alignment, args.max_dist)

    if args.moving_average_tsv:
        print_moving_average_tsv(seq1_aligned, seq2_aligned, args)
        return

    if rec1.description:
        print('>%s' % rec1.description)
    if rec2.description:
        print('>%s' % rec2.description)

    progress_1 = 0
    progress_2 = 0
    if args.start_pos:
        seq_num, start_pos = args.start_pos.split(":")
        seq_num = int(seq_num)
        start_pos = int(start_pos)
        if seq_num == 1:
            seq1_aligned, seq2_aligned, progress_1, progress_2 = trim_seqs(
                seq1_aligned, seq2_aligned, start_pos)
        else:
            seq2_aligned, seq1_aligned, progress_2, progress_1 = trim_seqs(
                seq2_aligned, seq1_aligned, start_pos)

    for seq1_line, seq2_line in zip(wrap(seq1_aligned, args.columns),
                                    wrap(seq2_aligned, args.columns)):
        if args.chart:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(10, 2), constrained_layout=True)

            xs = range(progress_1,
                       progress_1 + len(seq1_line))

            ax.plot(xs, list(seq1_line), label="Reference", color="blue")
            ax.plot(xs, list(seq2_line), label="Query", color="green")

            plt.show()
            plt.clf()

        seq1_line, seq2_line = color_mismatches(seq1_line, seq2_line)
        progress_1_start_str = "%s (%.0f%%)" % (
            progress_1, 100 * progress_1 / len(seq1))
        progress_1 += count_bases_in_printable_line(seq1_line)
        progress_1_end_str = "(%.0f%%) %s" % (
            100 * progress_1 / len(seq1), progress_1)
        if args.print_progress:
            print("%s%s%s" % (
                progress_1_start_str,
                " "*(max(1, len(rmcolor(seq1_line))
                     - len(progress_1_start_str)
                     - len(progress_1_end_str))),
                progress_1_end_str))

        print(seq1_line)
        print(seq2_line)
        progress_2_start_str = "%s (%.0f%%)" % (
            progress_2, 100 * progress_2 / len(seq2))
        progress_2 += count_bases_in_printable_line(seq2_line)
        progress_2_end_str = "(%.0f%%) %s" % (
            100 * progress_2 / len(seq2), progress_2)
        if args.print_progress:
            print("%s%s%s" % (
                progress_2_start_str,
                " "*(max(1, len(rmcolor(seq2_line))
                     - len(progress_2_start_str)
                     - len(progress_2_end_str))),
                progress_2_end_str))

        print()

def start():
    parser = argparse.ArgumentParser(
        description='Align sequences and show them vertically interleaved')
    parser.add_argument(
        'in1',
        metavar='[rc:]ACTG|path[:id]',
        help='First input.  Either a literal sequence or a path to a fasta '
        'or fastq file.  If a path, optionally include a colon-separated '
        'id to refer to identify a specific record in the file.')
    parser.add_argument(
        'in2', metavar='[:rc]ACTG|path[:id]', help='Same as first input.')
    parser.add_argument(
        '--columns', type=int, metavar='N',
        help='How many columns to wrap at.  If unspecified, autodetects.')
    parser.add_argument(
        '--max-dist', type=int, metavar='N', default=10,
        help='How hard to try to avoid over-alignment when blocks have been '
        'substituted out.')
    parser.add_argument(
        '--min-score', type=int, metavar='N', default=40,
        help='Minimum score of alignment to print.')
    parser.add_argument(
        '--print-progress', action='store_true',
        help='Mark lines with how far along the genomes they represent.')
    parser.add_argument(
        '--start-pos', type=str, metavar='[1:2]:N', default='',
        help='Start position, as bases into either sequence 1 or 1')
    parser.add_argument(
        '--moving-average-tsv', action='store_true',
        help='Instead of printing bases, output a tsv of how well these '
        'genomes align.')
    parser.add_argument(
        '--moving-average-width', type=int, metavar='N', default=40,
        help='Standard deviation of the moving average. '
        'Ignored unless --moving-average-tsv.')
    parser.add_argument(
        '--chart', action='store_true',
        help='visualize the alignment with some squiggles')
    args = parser.parse_args()

    if not args.columns:
        args.columns = get_columns()

    run(args)

if __name__ == '__main__':
    try:
        start()
    except BrokenPipeError:
        sys.stderr.close()
        pass
