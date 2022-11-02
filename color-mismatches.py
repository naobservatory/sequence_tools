import sys

from collections import Counter

COLOR_RED = '\x1b[1;31m'
COLOR_END = '\x1b[0m'

def start():
    out = [line.strip() for line in sys.stdin]
    max_out = max(len(x) for x in out)
    highlight = []
    for col in range(max_out):
        bases = Counter()
        for row in out:
            try:
                val = row[col]
            except IndexError:
                continue

            if val == ' ': continue
            bases[val] += 1

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

    for line in out:
        print(line)


if __name__ == "__main__":
    start()
