# fuzzy_highlighter

Highlight matching subsequences with fuzzy matching.  Handles substitutions but
not insertions or deletions.

## Installation

1. Install rust
2. `cargo build --release`
3. Put `target/release/fuzzy_highlighter` somewhere on your path.

## Usage

    cat input | fuzzy_highlighter pattern1 pattern2

Run `fuzzy_highlighter --help` for more details.
