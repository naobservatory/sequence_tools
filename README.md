# sequence tools

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

## fastq_viewer.py

```
usage: fastq_viewer.py [-h] [--colorize-quality] [--columns N] [--skip-lines]
                       [--highlighted-only] [--show-quality]
                       [--show-quality-when-highlighted] [--id-matches REGEX]
                       [--seq-matches REGEX[:COLOR]] [--max-quality CHAR]
                       [fastq_filenames ...]

Display FASTQ files in a more human-readable way.

positional arguments:
  fastq_filenames

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
  --seq-matches REGEX[:COLOR]
                        Only print sequences which match the regex. May be
                        specified multiple times, and sequences that match any
                        will be printed. Colors are red, yellow, green, blue,
                        magenta, cyan, white, and black, with an optional
                        _bold suffix.
  --max-quality CHAR    If your sequencer doesn't use the whole quality range,
                        set this to something smaller than '~' to make better
                        use of the available colors. Ignored when --colorize-
                        quality is false.
```
