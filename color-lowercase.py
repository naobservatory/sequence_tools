#!/usr/bin/env python3
import re
import sys

COLOR_RED = '\x1b[1;31m'
COLOR_END = '\x1b[0m'

for line in sys.stdin:
    if re.match("^[ACGTNagctn]+\n$", line):
        line = re.sub("([acgtn]+)", "%s\\1%s" % (COLOR_RED, COLOR_END), line)
    sys.stdout.write(line)
