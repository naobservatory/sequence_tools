#!/bin/bash
INVOCATION="./align.py"
PREFIX="align"
source test-util.sh

check align-all              tests/a.fasta tests/a.fasta
check align-specific         tests/a.fasta:seq123 tests/a.fasta:seq456
check align-fail             tests/a.fasta:seq123 tests/a.fasta:seq789
check min-score              tests/a.fasta:seq123 tests/a.fasta:seq789 \
                                --min-score 0
check align-fastq            tests/b.fastq tests/b.fastq
check literal                tests/b.fastq:fq2 ACAACCTACCGTGACAAAGAAAGTTGTC
