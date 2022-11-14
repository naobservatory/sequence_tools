#!/usr/bin/python3

# Usage: n2-align.py in out contig
#
# in format: reads or fasta
# out format: annotated fasta

import sys

MIN_OVERLAP=5

in_fname, out_fname, contig, = sys.argv[1:]

if len(contig) <= MIN_OVERLAP:
    raise Exception(
        "contig %s (length=%s) must be longer than MIN_OVERLAP (%s)" % (
            contig, len(contig), MIN_OVERLAP))

def rc(s):
    return "".join({'T':'A',
                    'G':'C',
                    'A':'T',
                    'C':'G',
                    'N':'N'}[x] for x in reversed(s))

def score_at(contig, read, contig_pos):
    score = 0
    for i in range(len(read)):
        read_val = read[i]
        contig_val = None
        if 0 <= contig_pos + i < len(contig):
            contig_val = contig[contig_pos + i]
        if read_val == contig_val:
            score += 1
    return score

def align_read(read):
    first_eligible_position = -(len(read) - MIN_OVERLAP)
    last_eligible_position = len(contig) - MIN_OVERLAP
    best_score = 0
    best_contig_pos = None
    for contig_pos in range(first_eligible_position,
                            last_eligible_position+1):
        score = score_at(contig, read, contig_pos)
        if score > best_score:
            best_score = score
            best_contig_pos = contig_pos
        
    return best_contig_pos, best_score

with open(in_fname) as inf:
    with open(out_fname, "w") as outf:
        seqid = ""
        for line in inf:
            line = line.strip()
            if line.startswith(">"):
                seqid = line[1:]
                continue
            
            pos, score = align_read(line)
            rc_line = rc(line)
            rc_pos, rc_score = align_read(rc_line)

            if rc_score > score:
                line = rc_line
                pos = rc_pos

            outf.write(">%s%spos=%s\n%s\n" % (
                seqid,
                " "if seqid else "",
                pos,
                line))
            seqid = ""
