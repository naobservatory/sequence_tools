# sequence tools

## align.py

```
usage: align.py [-h] [--columns N] [--max-dist N] [--min-score N]
                ACTG|path[:id] in2

Align sequences and show them vertically interleaved

positional arguments:
  ACTG|path[:id]  First input. Either a literal sequence or a path to a fasta
                  or fastq file. If a path, optionally include a colon-
                  separated id to refer to identify a specific record in the
                  file.
  in2             See in1

optional arguments:
  -h, --help      show this help message and exit
  --columns N     How many columns to wrap at. If unspecified, autodetects.
  --max-dist N    How hard to try to avoid over-alignment when blocks have
                  been substituted out.
  --min-score N   Minimum score of alignment to print.
```

## seqdsp.py

```
usage: seqdsp.py [-h] [--colorize-quality] [--columns N] [--skip-lines]
                 [--highlighted-only] [--show-quality]
                 [--show-quality-when-highlighted] [--id-matches REGEX]
                 [--seq-matches [rc]ACTG[:COLOR]] [--max-quality CHAR]
                 [--min-score N]
                 [fname[:id] ...]

Display seqence files in a human-readable way.

positional arguments:
  fname[:id]

optional arguments:
  -h, --help            show this help message and exit
  --colorize-quality    Color the quality line to show the quality visually.
                        Order, lowest to highest, is black, magenta, red,
                        yellow, white, green, blue, cyan.
  --columns N           How many columns to wrap at. If unspecified,
                        autodetects.
  --skip-lines          Leave extra space between lines for readability.
  --highlighted-only    Only print sequences that have colored portions. For
                        example, the output of fuzzy_highlighter.
  --show-quality        Print quality lines.
  --show-quality-when-highlighted
                        Print quality lines corresponding to sequence lines
                        that have colored portions.
  --id-matches REGEX    Only print sequences whose id line matches the regex.
  --seq-matches [rc]ACTG[:COLOR]
                        Only print matching sequence, or reverse complement if
                        "rc" prefix is present. May bespecified multiple
                        times, and sequences that match any will be printed.
                        Colors are red, yellow, green, blue, magenta, cyan,
                        white, and black.
  --max-quality CHAR    If your sequencer doesn't use the whole quality range,
                        set this to something smaller than '~' to make better
                        use of the available colors. Ignored when --colorize-
                        quality is false.
  --min-score N         Minimum score of alignment to print.
```

## fuzzy_highlighter

Highlight matching subsequences with fuzzy matching.  Handles substitutions but
not insertions or deletions.

### Installation

1. Install rust
2. `cargo build --release`
3. Put `target/release/fuzzy_highlighter` somewhere on your path.

### Usage

    cat input | fuzzy_highlighter pattern1 pattern2

Run `fuzzy_highlighter --help` for more details.

