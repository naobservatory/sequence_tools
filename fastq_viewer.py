#!/usr/bin/env python3

import argparse
import fileinput
import sys
import re
import os

from Bio import Align
from Bio.Seq import Seq

import ansiwrap # python3 -m pip install ansiwrap

from util import COLORS, COLOR_END
from util import die, guess_format_or_die, get_columns
from util import interpret_sequence_argument

def colorless_len(line):
    l1 = ansiwrap.ansilen(line)
    l2 = len(ansistrip(line))
    if l1 != l2:
        die('Unable to determine length.  Got both %s and %s for "%s"' % (
            l1, l2, line))
    return l1

def ansistrip(line):
    # This is much more agressive than ANSIRE, but it should be the same for
    # our use case.  The key difference is at end of line with an incomplete
    # escape code when ANSIRE won't strip and we will.

    # Nucleotides are generally in ACGT but other letters are used
    # occasionally and we don't need to be picky.
    return re.sub('[^A-Z]+', '', line)

def has_color(line):
    return re.search('[^A-Z]+', line)

COLOR_BUCKETS = [
    COLORS['black_bold'],  # bold so it's not invisible
    COLORS['magenta'],
    COLORS['red'],
    COLORS['yellow'],
    COLORS['white'],
    COLORS['green'],
    COLORS['blue'],
    COLORS['cyan'],
]

def colorize_quality(line, max_quality):
    # Quality is ascii 33 (!) through 126 (~) in order, very tidy.
    out = []
    for c in line:
        if c < '!' or c > '~':
            raise Exception('Value %r out of range in %r' % (
                c, line))
        bucket = int(len(COLOR_BUCKETS) *
                     (ord(c) - ord('!')) / (ord(max_quality) - ord('!')))
        bucket = min(len(COLOR_BUCKETS)-1, bucket)
        out.append('%s%s%s' % (
            COLOR_BUCKETS[bucket], c, COLOR_END))

    return ''.join(out)

def wrap(text, strip_ansi, columns):
    # You would think we could just use ansiwrap.wrap, but that doesn't
    # actually work: https://github.com/jonathaneunice/ansiwrap/issues/10
    # Do it ourselves.
    #
    # This is easier than the real way, because nucleotide sequences are always
    # [A-Z]+ and escape sequences are fully out of band.
    line = ''
    lines = []
    for c in text:
        line += c

        length = len(ansistrip(line)) if strip_ansi else len(line)
        if length == columns:
            lines.append(line)
            line = ''
    if line:
        lines.append(line)

    if strip_ansi:
        return ansiwrap.ansi_terminate_lines(lines)
    else:
        return lines

def rc(s):
    return "".join({'T':'A',
                    'G':'C',
                    'A':'T',
                    'C':'G'}[x] for x in reversed(s))

def print_record(record, args):
    if args.id_matches and not re.search(args.id_matches, record.description):
        return

    seq = str(record.seq)

    if args.seq_matches:
        aligner = Align.PairwiseAligner()
        # These are the scoring settings porechop uses by default.
        # https://github.com/rrwick/Porechop/blob/master/porechop/porechop.py#L145
        aligner.end_gap_score = 0
        aligner.match_score = 3
        aligner.mismatch_score = -6
        aligner.internal_open_gap_score = -5
        aligner.internal_extend_gap_score = -2

        any_matched = False
        for seq_matcher in args.seq_matches:
            color = 'red'
            if ':' in seq_matcher:
                seq_matcher, color = seq_matcher.split(':')

            if seq_matcher.startswith('rc'):
                seq_matcher = rc(seq_matcher[2:])

            if color:
                if color not in COLORS:
                    die('Unknown color %r' % color)
                else:
                    color = color + "_bold"

            if aligner.score(seq, seq_matcher) / min(
                    len(seq), len(seq_matcher)) < (args.min_score/100):
                continue

            alignment_seq, _ = aligner.align(seq, seq_matcher)[0].aligned
            new_seq = []
            pos = 0
            for start, end in alignment_seq:
                new_seq.append(seq[pos:start])
                new_seq.append(COLORS[color])
                new_seq.append(seq[start:end])
                new_seq.append(COLOR_END)
                pos = end
            new_seq.append(seq[pos:])
            seq = ''.join(new_seq)
            any_matched = True

        if not any_matched:
            return

    if args.highlighted_only and not has_color(seq):
        return

    print(COLORS['black_bold'] + record.description + COLOR_END)

    if 'phred_quality' in record.letter_annotations:
        quality =''.join(chr(ord('!') + x)
                         for x in record.letter_annotations['phred_quality'])
        for sequence_line, quality_line in zip(
                wrap(seq, strip_ansi=True, columns=args.columns),
                wrap(quality, strip_ansi=False, columns=args.columns)):
            if ansiwrap.ansilen(sequence_line) != ansiwrap.ansilen(quality_line):
                print(sequence_line)
                print(quality_line)
                raise Exception('%s vs %s' % (
                    ansiwrap.ansilen(sequence_line),
                    ansiwrap.ansilen(quality_line)))

            print(sequence_line)
            if args.show_quality or (args.show_quality_when_highlighted and
                                     has_color(sequence_line)):
                if args.colorize_quality:
                    quality_line = colorize_quality(quality_line,
                                                    max_quality=args.max_quality)
                print(quality_line)
                if args.skip_lines:
                    print()
    else:
        print(seq)

def start():
    parser = argparse.ArgumentParser(
        description='Display FASTQ files in a more human-readable way.')
    parser.add_argument('fastq_filenames', nargs='*', metavar='fname[:id]')
    parser.add_argument(
        '--colorize-quality', action='store_true',
        help='Color the quality line to show the quality visually.  Order, '
        'lowest to highest, is black, magenta, red, yellow, white, green, '
        'blue, cyan.')
    parser.add_argument(
        '--columns', type=int, metavar='N',
        help='How many columns to wrap at.  If unspecified, autodetects.')
    parser.add_argument(
        '--skip-lines', action='store_true',
        help='Leave extra space between lines for readability.')
    parser.add_argument(
        '--highlighted-only', action='store_true',
        help='Only print sequences that have colored portions. For example, '
        'the output of fuzzy_highlighter.')
    parser.add_argument(
        '--show-quality', action='store_true', help='Print quality lines.')
    parser.add_argument(
        '--show-quality-when-highlighted', action='store_true',
        help='Print quality lines corresponding to sequence lines that '
        'have colored portions.')
    parser.add_argument(
        '--id-matches', metavar='REGEX',
        help='Only print sequences whose id line matches the regex.')
    parser.add_argument(
        '--seq-matches', metavar='[rc]ACTG[:COLOR]',
        action='append',
        help='Only print matching sequence, or reverse complement if "rc" '
        'prefix is present. May bespecified multiple times, and sequences '
        'that match any will be printed. Colors are red, yellow, green, blue, '
        'magenta, cyan, white, and black.')
    parser.add_argument(
        '--max-quality', metavar='CHAR', default='D',
        help="If your sequencer doesn't use the whole quality range, set this "
        "to something smaller than '~' to make better use of the available "
        "colors.  Ignored when --colorize-quality is false.")
    parser.add_argument(
        '--min-score', type=int, metavar='N', default=40,
        help='Minimum score of alignment to print.')

    args = parser.parse_args()

    if (len(args.max_quality) != 1 or
        args.max_quality < '!' or
        args.max_quality > '~'):
        die('--max-quality must be between "!" and "~"; got "%s"' %
            args.max_quality)

    if not args.columns:
        args.columns = get_columns()

    run(args)

def run(args):
    for fname in args.fastq_filenames:
        for record in interpret_sequence_argument(fname):
            print_record(record, args)

if __name__ == '__main__':
    start()
