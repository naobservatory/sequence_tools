from Bio import SeqIO
from Bio.Seq import Seq
from collections import namedtuple
from collections import defaultdict
import regex
import glob
import itertools
import json
import sys

Barcode = namedtuple('Barcode', ['barcode_id', 'guppy_barcode_id', 'fwd', 'rev'])

def parse_barcodes():
    barcodes = []
    with open("bararr.csv") as inf:
        for i, line in enumerate(inf):
            line = line.strip()
            cols = line.split("\t")
            if len(cols) == 4:
                barcode_id, _, fwd, rev = cols
                barcodes.append(Barcode(
                    barcode_id, 'bw_%s' % (str(len(barcodes)).zfill(3)),
                    Seq(fwd), Seq(rev)))
    return barcodes

def contains(haystack, needle, max_errors):
    if not max_errors:
        return needle in haystack
    return regex.search('(%s){e<=%s}' % (needle, max_errors), str(haystack))

def demultiplex(record, barcodes, max_errors):
    merr = max_errors
    matches = []
    for barcode in barcodes:
        fwd = barcode.fwd
        rev = barcode.rev
        fwd_rc = fwd.reverse_complement()
        rev_rc = rev.reverse_complement()

        head = record.seq[:100]
        tail = record.seq[-100:]

        if contains(head, fwd, merr) and contains(tail, rev_rc, merr):
            matches.append(barcode.guppy_barcode_id)
        elif contains(head, rev, merr) and contains(tail, fwd_rc, merr):
            matches.append(barcode.guppy_barcode_id)

    return matches[0] if len(matches) == 1 else "unclassified"

def start(fnames):
    barcodes = parse_barcodes()

    for fname in fnames:
        with open(fname) as inf:
            for record in SeqIO.parse(fname, "fastq"):
                matched_6 = False
                for max_errors in [6,7]:
                    our_barcode = demultiplex(record, barcodes, max_errors)
                    matched = our_barcode != "unclassified"

                    if matched:
                        if max_errors == 6:
                            matched_6 = True
                        elif not matched_6:
                            print ("%s:%s\t%s" % (fname, record.id, our_barcode))

if __name__ == "__main__":
    start(sys.argv[1:])
