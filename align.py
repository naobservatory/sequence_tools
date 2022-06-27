import os
import argparse

from Bio import Align

"""
aligner = Align.PairwiseAligner()
alignment = aligner.align(a, b)[0]

c = str(alignment).split('\n')[0]
d = str(alignment).split('\n')[2]

c = c.replace('-', ' ')
d = d.replace('-', ' ')

COLOR_RED = '\x1b[1;31m'
COLOR_END = '\x1b[0m'

for e, f in zip(wrap(c), wrap(d)):
    g = ''
    h = ''
    for w, x in zip(e.ljust(35), f.ljust(35)):
        if w == x:
            g += x
            h += x
        else:
            g += COLOR_RED + x + COLOR_END
            h += COLOR_RED + w + COLOR_END
    
    print ("%s   %s" % (g, h))
"""

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

def run(args):
    aligner = Align.PairwiseAligner()
    alignment = aligner.align(args.seq1, args.seq2)[0]
    seq1_aligned, _, seq2_aligned, _ = str(alignment).split('\n')

    seq1_aligned, seq2_aligned = collapse_subs(
        seq1_aligned, seq2_aligned, args.max_dist)
    
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
    parser.add_argument('seq1')
    parser.add_argument('seq2')
    parser.add_argument(
        '--columns', type=int, metavar='N',
        help='How many columns to wrap at.  If unspecified, autodetects.')
    parser.add_argument(
        '--max-dist', type=int, metavar='N', default=10,
        help='How hard to try to avoid over-alignment when blocks have been '
        'substituted out.')
    args = parser.parse_args()
    
    if not args.columns:
        args.columns = get_columns()

    run(args)
    
if __name__ == '__main__':
    start()
