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
                    barcode_id, 'bw_%s' % (str(len(barcodes)).zfill(3)),
                    Seq(fwd), Seq(rev)))
    return barcodes

def contains(haystack, needle):
    return regex.search('(%s){e<=6}' % needle, str(haystack))


def demultiplex(record, barcodes):
    matches = []
    for barcode in barcodes:
        fwd = barcode.fwd
        rev = barcode.rev
        fwd_rc = fwd.reverse_complement()
        rev_rc = rev.reverse_complement()

        head = record.seq[:100]
        tail = record.seq[-100:]

        if contains(head, fwd) and contains(tail, rev_rc):
            matches.append(barcode.guppy_barcode_id)
        elif contains(head, rev) and contains(tail, fwd_rc):
            matches.append(barcode.guppy_barcode_id)

    return matches[0] if len(matches) == 1 else "unclassified"
            
def start():
    barcodes = parse_barcodes()

    n_success = 0
    n_total = 0
    for fname in glob.glob("pass/**/*.fastq"):
        n_file_success = 0
        n_file_total = 0
        with open(fname) as inf:
            for record in SeqIO.parse(fname, "fastq"):
                our_barcode = demultiplex(record, barcodes)
                guppy_barcode = record.description.split()[-1].removeprefix(
                    "barcode=")
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
