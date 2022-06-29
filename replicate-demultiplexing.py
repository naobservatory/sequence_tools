from Bio import SeqIO
from Bio.Seq import Seq
from collections import namedtuple
import regex
import glob

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
                    barcode_id, 'bw_%s' % (str(len(barcodes)-1).zfill(3)),
                    Seq(fwd), Seq(rev)))
    return barcodes

def contains(haystack, needle):
    max_errors = len(needle) // 4
    #return regex.search('(%s){e<=%s}' % (needle, max_errors), str(haystack))
    return regex.search('(%s){e<=4}' % (needle), str(haystack))

max_distance_from_end = 20
def head_contains(seq, needle):
    return contains(seq[:len(needle)+max_distance_from_end], needle)

def tail_contains(seq, needle):
    return contains(seq[-len(needle)+max_distance_from_end:],
                    needle.reverse_complement())

def demultiplex(record, barcodes):
    matches = []
    for barcode in barcodes:
        if ((head_contains(record.seq, barcode.fwd) and
             tail_contains(record.seq, barcode.rev)) or
            (head_contains(record.seq, barcode.rev) and
             tail_contains(record.seq, barcode.fwd))):
            matches.append(barcode.guppy_barcode_id)

    return matches[0] if len(matches) == 1 else "unclassified"
            
def start():
    barcodes = parse_barcodes()

    n_success = 0
    n_total = 0
    #for fname in glob.glob("pass/**/*.fastq"):
    for fname in ['pass/bw_063/fastq_runid_1b7bf189bb6e7b6627c717b20111ff1777ce9409_22_0.fastq']:
        n_file_success = 0
        n_file_total = 0
        with open(fname) as inf:
            for record in SeqIO.parse(fname, "fastq"):
                our_barcode = demultiplex(record, barcodes)
                guppy_barcode = record.description.split()[-1].removeprefix(
                    "barcode=")
                print("%s %s %s" % (record.id, our_barcode, guppy_barcode))
                success = our_barcode == guppy_barcode
                if success:
                    n_success += 1
                    n_file_success += 1

                n_total += 1
                n_file_total += 1
        print("%.0f%%: %s" % (
            n_file_success/n_file_total*100, fname))
    print("%.0f%%: %s / %s" % (n_success/n_total*100, n_success, n_total))

if __name__ == "__main__":
    start()
