#!/usr/bin/env python3
import re
import sys
import argparse
from collections import Counter

COLOR_RED = '\x1b[1;31m'
COLOR_CYAN = '\x1b[1;36m'
COLOR_END = '\x1b[0m'

def remove_count_prefixes_if_present(lines):
    counts = []
    out = []
    for line in lines:
        if line.count("\t") != 1:
            break
        count, line = line.split("\t")
        if not count.isdigit():
            break
        out.append(line)
        counts.append(int(count))
    else:
        # completed successfully
        return counts, out

    # not the right format, just put every count as equally common
    return [1]*len(lines), lines

def start():
    parser = argparse.ArgumentParser(
        description='Given horizontally aligned sequences,'
                    ' highlight disagreement')
    parser.add_argument(
        '--min-count', type=int, default=0, metavar='N',
        help="Sequences appearing less often won't be printed")
    run(parser.parse_args())

def run(args):
    # remove final linebreak
    out = [line[:-1] for line in sys.stdin]

    counts, out = remove_count_prefixes_if_present(out)

    if out and re.match("^[0-9 ]*$", out[0]):
        print(COLOR_CYAN + out[0] + COLOR_END)
        out = out[1:]
    
    max_out = max(len(x) for x in out)
    highlight = []
    for col in range(max_out):
        bases = Counter()
        for count, row in zip(counts, out):
            try:
                val = row[col]
            except IndexError:
                continue

            if val == ' ': continue
            bases[val] += count

        most_common_val = None
        most_common_vals = bases.most_common(1)
        if most_common_vals:
            (most_common_val, most_common_count), = most_common_vals

            for base, count in bases.items():
                if base != most_common_val and count == most_common_count:
                    most_common_val = None
                    break

        for row, line in enumerate(out):
            try:
                val = line[col]
            except IndexError:
                continue

            if val == ' ': continue

            if val != most_common_val:
                highlight.append((row, col))

    for row, col in reversed(highlight):
        out[row] = "%s%s%s%s%s" % (
            out[row][:col],
            COLOR_RED,
            out[row][col],
            COLOR_END,
            out[row][col+1:])

    for count, line in zip(counts, out):
        if count < args.min_count: continue
        print(line)


if __name__ == "__main__":
    start()
