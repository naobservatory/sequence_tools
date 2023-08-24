import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord

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
BEGIN_UNDERLINE = '\x1b[4m'
END_UNDERLINE = '\x1b[24m'

def die(msg):
    print(COLORS['red'] + msg + COLOR_END)
    sys.exit(1)

def guess_format_or_die(fname):
    _, ext = os.path.splitext(fname)
    if ext in ['.fasta', '.fa', '.faa']:
        return 'fasta'
    elif ext in ['.fastq', '.fq']:
        return 'fastq'
    else:
        die('Unable to guess format of %r from extension' % fname)

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

def interpret_sequence_argument(arg_in):
    is_rc = False
    if arg_in.startswith("rc:"):
        arg_in = arg_in.removeprefix("rc:")
        is_rc = True

    if ':' in arg_in:
        seq, seq_id = arg_in.split(':')
    else:
        seq = arg_in
        seq_id = None

    def maybe_rc(record):
        if is_rc:
            record.seq = record.seq.reverse_complement()
        return record

    if os.path.exists(seq):
        # sequences on disk
        fname = seq
        with open(fname) as inf:
            records = SeqIO.parse(inf, guess_format_or_die(fname))
            if seq_id:
                for record in records:
                    if record.id == seq_id:
                        return [maybe_rc(record)]
                else:
                    die('Sequence %r not found in %r' % (
                        seq_id, fnamer))
            else:
                return [maybe_rc(record) for record in records]
    else:
        # raw on the command line
        return [maybe_rc(SeqRecord(Seq(seq), description=''))]
