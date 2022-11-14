#!/usr/bin/python3

# Usage: cat reads | display-reads.py contig

import sys
import re

MIN_OVERLAP=10

contig, = sys.argv[1:]

if len(contig) <= MIN_OVERLAP:
    raise Exception(
        "contig (length=%s) must be longer than MIN_OVERLAP (%s)" % (
            len(contig), MIN_OVERLAP))

reads = []
for line in sys.stdin:
    line = line.strip()
    if line.startswith(">"):
        seqid = line[1:]
        continue

    pos, = re.findall("pos=([-0-9]+)", seqid)
    pos = int(pos)
    
    reads.append((int(pos), line))
    seqid = None

if not reads:
    raise Exception("no reads -- supply on stdin")
    
reads.sort()

contig_offset = -reads[0][0]
if contig_offset < 0:
    contig_offset = 0

print("%s%s" % (
    " "*contig_offset, contig))

for contig_pos, read in reads:
    print("%s%s" % (
        " "*(contig_pos + contig_offset), read))

