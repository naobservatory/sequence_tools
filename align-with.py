import sys

from collections import Counter

COLOR_RED = '\x1b[1;31m'
COLOR_END = '\x1b[0m'


pos_lines = []
out = []
def start(needle):
    for line in sys.stdin:
        if needle not in line: continue

        pos = line.index(needle)

        pos_lines.append((pos, line.strip()))

    if not pos_lines: return

    pos_lines.sort(reverse=True)
    max_pos = pos_lines[0][0]

    for pos, line in pos_lines:
        out.append("%s%s" % (" "*(max_pos-pos), line))

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
    start(*sys.argv[1:])
