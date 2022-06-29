#!/usr/bin/env python3

import os
import argparse
import itertools
import sys

from Bio import Align
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord

COLOR_RED = '\x1b[1;31m'
COLOR_GREEN = '\x1b[1;32m'
COLOR_YELLOW = '\x1b[1;33m'
COLOR_END = '\x1b[0m'

def die(msg):
    print(msg)
    sys.exit(1)

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

# copied from icdiff
def get_columns():
    if os.name == 'nt':
        try:
            import struct
            from ctypes import windll, create_string_buffer

            fh = windll.kernel32.GetStdHandle(-12)  # stderr is -12
            csbi = create_string_buffer(22)
            windll.kernel32.GetConsoleScreenBufferInfo(fh, csbi)
            res = struct.unpack('hhhhHhhhhhh', csbi.raw)
            return res[7] - res[5] + 1  # right - left + 1

        except Exception:
            pass

    else:

        def ioctl_GWINSZ(fd):
            try:
                import fcntl
                import termios
                import struct

                cr = struct.unpack(
                    'hh', fcntl.ioctl(fd, termios.TIOCGWINSZ, '1234')
                )
            except Exception:
                return None
            return cr

        cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
        if cr and cr[1] > 0:
            return cr[1]
    return 80

def collapse_subs(seq1, seq2, max_dist):
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

def color_mismatches(seq1_line, seq2_line):
    out1 = []
    out2 = []

    for c1, c2 in zip(seq1_line, seq2_line):
        c1_color = None
        c2_color = None
        if c1 == '-':
            c2_color = COLOR_GREEN
        elif c2 == '-':
            c1_color = COLOR_RED
        elif c1 != c2:
            c1_color = COLOR_YELLOW
            c2_color = COLOR_YELLOW

        if c1_color:
            out1.append(c1_color)
        out1.append(c1)
        if c1_color:
            out1.append(COLOR_END)

        if c2_color:
            out2.append(c2_color)
        out2.append(c2)
        if c2_color:
            out2.append(COLOR_END)

    return ''.join(out1), ''.join(out2)

def interpret_sequence_argument(arg_in):
    if ':' in arg_in:
        seq, seq_id = arg_in.split(':')
    else:
        seq = arg_in
        seq_id = None

    if os.path.exists(seq):
        # sequences on disk
        with open(seq) as inf:
            if seq.endswith('.fastq'):
                records = SeqIO.parse(inf, 'fastq')
            elif seq.endswith('.fasta') or seq.endswith('.fa'):
                records = SeqIO.parse(inf, 'fasta')
            else:
                die('unknown file format for %r' % seq)

            if seq_id:
                for record in records:
                    if record.id == seq_id:
                        return [record]
                else:
                    die('Sequence %r not found in %r' % (
                        seq_id, seq))
            else:
                return list(records)
    else:
        # raw on the command line
        return [SeqRecord(Seq(seq), description='')]

def run(args):
    recs1 = interpret_sequence_argument(args.in1)
    recs2 = interpret_sequence_argument(args.in2)

    for rec1, rec2 in itertools.product(recs1, recs2):
        align_and_print(rec1, rec2, args)

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

    seq1_aligned, _, seq2_aligned, _ = str(alignment).split('\n')

    seq1_aligned, seq2_aligned = collapse_subs(
        seq1_aligned, seq2_aligned, args.max_dist)

    if rec1.description:
        print('>%s' % rec1.description)
    if rec2.description:
        print('>%s' % rec2.description)

    for seq1_line, seq2_line in zip(
            wrap(seq1_aligned, args.columns),
            wrap(seq2_aligned, args.columns)):
        seq1_line, seq2_line = color_mismatches(seq1_line, seq2_line)
        print(seq1_line)
        print(seq2_line)
        print()

def start():
    parser = argparse.ArgumentParser(
        description='Align sequences and show them vertically interleaved')
    parser.add_argument(
        'in1',
        metavar='ACTG|path[:id]',
        help='First input.  Either a literal sequence or a path to a fasta '
        'or fastq file.  If a path, optionally include a colon-separated '
        'id to refer to identify a specific record in the file.')
    parser.add_argument('in2', help='See in1')
    parser.add_argument(
        '--columns', type=int, metavar='N',
        help='How many columns to wrap at.  If unspecified, autodetects.')
    parser.add_argument(
        '--max-dist', type=int, metavar='N', default=10,
        help='How hard to try to avoid over-alignment when blocks have been '
        'substituted out.')
    parser.add_argument(
        '--min-score', type=int, metavar='N', default=40,
        help='Minimum score of alignment to print, counting matches as +1, '
        'mismatches as -1, and gaps as -1.')
    args = parser.parse_args()

    if not args.columns:
        args.columns = get_columns()

    run(args)

if __name__ == '__main__':
    start()
