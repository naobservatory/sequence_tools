import json

MAX_ERRORS=8

with open("demultiplexed-by-error-allowance.jsons") as inf:
    results = [{
        "totals": 0,
        "classified": 0,
        "unclassified": 0,
    } for i in range(MAX_ERRORS)]
    for line in inf:
        fname, counts = json.loads(line)
        for max_edits, barcode_counts in enumerate(counts):
            for barcode, count in barcode_counts[0].items():
                results[max_edits]["totals"] += count
                if barcode == "unclassified":
                    results[max_edits]["unclassified"] += count
                else:
                    results[max_edits]["classified"] += count
    print("max_edits\tclassified\ttotal")
    for max_edits, result in enumerate(results):
        print("%s\t%s\t%s" % (
            max_edits,
            result["classified"],
            result["totals"]))
