#!/usr/bin/env python3

"""
Usage: ./kraken-highlight.py ACTG... 0:144 2697049:5 0:7 2045208:4 0:89
"""


import sys

K=35
seq, *krakens = sys.argv[1:]

taxids = set()

# start, length, taxid
parsed_krakens = []
start = 0
for kraken in krakens:
    taxid, length = kraken.split(":")
    length = int(length)
    parsed_krakens.append((start, length, taxid))
    taxids.add(taxid)
    start += length

if len(seq) != start+K-1:
    print("length mismatch: sequence is %s bases but with K=%s kraken "
          "output only covers %s" % (len(seq), K, start+K-1))

max_taxid_len = max(len(x) for x in taxids)

print("%s %s" % (
    "seq".rjust(max_taxid_len),
    seq))
for start, length, taxid in parsed_krakens:
    if taxid == "0":
        continue
    print("%s %s%s" % (
        taxid.rjust(max_taxid_len),
        " "*start,
        seq[start:start+length+K]))
      
