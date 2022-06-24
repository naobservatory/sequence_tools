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
$ python3 fastq_viewer.py  --help
usage: fastq_viewer.py [-h] [--colorize-quality] [--columns COLUMNS]
                       [--skip-lines] [--max-quality MAX_QUALITY]
                       [fastq_filenames ...]

Display FASTQ files in a more human-readable way. Wraps lines to terminal
width, preserving existing ansi colors, and interleaves sequence and quality
lines.

positional arguments:
  fastq_filenames

optional arguments:
  -h, --help            show this help message and exit
  --colorize-quality    Color the quality line to show the quality visually.
                        Order, lowest to highest, is black, magenta, red,
                        yellow, white, green, blue, cyan.
  --columns COLUMNS     How many columns to wrap at. If unspecified we
                        autodetect.
  --skip-lines          Leave extra space between lines for readability.
  --max-quality MAX_QUALITY
                        If your sequencer doesn't use the whole quality range,
                        set this to something smaller than '~' to make better
                        use of the available colors. Ignored when
                        --colorize-quality is false.
```
