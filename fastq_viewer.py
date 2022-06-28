import argparse
import fileinput
import sys
import re
import os

import ansiwrap # python3 -m pip install ansiwrap
import regex # python3 -m pip install regex

def die(msg):
    print(msg)
    sys.exit(1)

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

COLORS = {
    'black': '\x1b[30m',
    'red': '\x1b[31m',
    'green': '\x1b[32m',
    'yellow': '\x1b[33m',
    'blue': '\x1b[34m',
    'magenta': '\x1b[35m',
    'cyan': '\x1b[36m',
    'white': '\x1b[37m',

    'black_bold': '\x1b[1;30m',
    'red_bold': '\x1b[1;31m',
    'green_bold': '\x1b[1;32m',
    'yellow_bold': '\x1b[1;33m',
    'blue_bold': '\x1b[1;34m',
    'magenta_bold': '\x1b[1;35m',
    'cyan_bold': '\x1b[1;36m',
    'white_bold': '\x1b[1;37m',
}
COLOR_END = '\x1b[0m'


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

# copied from icdiff
def get_columns():
    if os.name == "nt":
        try:
            import struct
            from ctypes import windll, create_string_buffer

            fh = windll.kernel32.GetStdHandle(-12)  # stderr is -12
            csbi = create_string_buffer(22)
            windll.kernel32.GetConsoleScreenBufferInfo(fh, csbi)
            res = struct.unpack("hhhhHhhhhhh", csbi.raw)
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
                    "hh", fcntl.ioctl(fd, termios.TIOCGWINSZ, "1234")
                )
            except Exception:
                return None
            return cr

        cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
        if cr and cr[1] > 0:
            return cr[1]
    return 80

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

def print_seq(args, id_line, sequence, plus_line=None, quality=None):
    if args.highlighted_only and not has_color(sequence):
        return

    if args.id_matches and not regex.search(args.id_matches, id_line):
        return

    if args.seq_matches:
        any_matched = False
        for seq_matcher in args.seq_matches:
            color = None
            if seq_matcher.count(':') == 2:
                literal_seq, max_errors, color = seq_matcher.split(':')
                if not max_errors.isdigit():
                    die('Bad max_errors %r in %r' % (
                        max_errors, seq_matcher))
                seq_matcher = '(%s){e<=%s}' % (literal_seq, max_errors)
            elif seq_matcher.count(':') == 1:
                seq_matcher, color = seq_matcher.split(':')
            elif ':' in seq_matcher:
                die('Too many colons in %r' % seq_matcher)

            if color:
                if color not in COLORS:
                    die('Unknown color %r' % color)
                else:
                    color = color + "_bold"

            s = regex.search(seq_matcher, sequence)
            if s:
                any_matched = True
                start, end = s.span()
                if color:
                    sequence = (
                        sequence[:start] + COLORS[color] +
                        sequence[start:end] + COLOR_END + sequence[end:])
        if not any_matched:
            return

    print(COLORS['black_bold'] + id_line + COLOR_END)

    if quality:
        for sequence_line, quality_line in zip(
                wrap(sequence, strip_ansi=True, columns=args.columns),
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
        print(sequence)
        #for sequence_line in wrap(sequence, strip_ansi=True,
        #                          columns=args.columns):
        #    print(sequence_line)

    if plus_line:
        print(plus_line)

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
        '--seq-matches', metavar='(REGEX|ACTG:EDITS)[:COLOR]', action='append',
        help='Only print sequences which match the regex or sequence.  May be '
        'specified multiple times, and sequences that match any will be '
        'printed. Colors are red, yellow, green, blue, magenta, cyan, white, '
        'and black.')
    parser.add_argument(
        '--max-quality', metavar='CHAR', default='D',
        help="If your sequencer doesn't use the whole quality range, set this "
        "to something smaller than '~' to make better use of the available "
        "colors.  Ignored when --colorize-quality is false.")
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
        seq_id = None
        if ':' in fname:
            fname, seq_id = fname.split(':')
        with open(fname) as inf:
            view_file(args, inf, seq_id)

def view_file(args, inf, seq_id):
    format = None  # '@' for fastq or '>' for fasta

    id_line = None
    plus_line = None
    sequence = []
    quality = []
    quality_len = 0
    sequence_len = 0

    for lineno, line in enumerate(inf):
        line = line.strip()

        if not format:
            if line.startswith('>'):
                format = 'fasta'
            elif line.startswith('@'):
                format = 'fastq'
            else:
                die('%s:%s: unknown format %r' %
                    fileinput.filename(), lineno, line)

        if format == 'fastq':
            if not id_line:
                if not line.startswith('@'):
                    die('%s:%s: bad value: %r' % (
                        fileinput.filename(), lineno, line))
                id_line = line
            elif not plus_line:
                if line.startswith('+'):
                    plus_line = line
                else:
                    sequence.append(line)
                    sequence_len += colorless_len(line)
            else:
                quality.append(line)
                quality_len += len(line)

                if quality_len == sequence_len:
                    if (not seq_id or
                        id_line[1:] == seq_id  or
                        id_line[1:].startswith(seq_id + " ")):

                        print_seq(args,
                                  id_line,
                                  ''.join(sequence),
                                  plus_line,
                                  quality=''.join(quality))
                    id_line = None
                    plus_line = None
                    sequence = []
                    quality = []
                    quality_len = 0
                    sequence_len = 0

                elif quality_len > sequence_len:
                    die('quality longer than sequence at line %s of %s '
                        '(%s > %s)' % (lineno, fileinput.filename(), quality_len,
                                       sequence_len))
        elif format == 'fasta':
            if line.startswith('>'):
                if sequence:
                    print_seq(args, id_line, ''.join(sequence))
                    id_line = None
                    sequence = []
                id_line = line
            else:
                sequence.append(line)
    if format == 'fasta' and sequence:
        print_seq(args, id_line, ''.join(sequence))



if __name__ == '__main__':
    start()
