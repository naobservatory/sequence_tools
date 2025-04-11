#!/bin/bash
INVOCATION="./align.py --columns=80"
PREFIX="align"
source test-util.sh

check align-all              tests/a.fasta tests/a.fasta
check align-specific         tests/a.fasta:seq123 tests/a.fasta:seq456
check align-specific-num     tests/a.fasta:0 tests/a.fasta:1
check align-fail             tests/a.fasta:seq123 tests/a.fasta:seq789
check min-score              tests/a.fasta:seq123 tests/a.fasta:seq789 \
                                --min-score 0
check align-fastq            tests/b.fastq tests/b.fastq
check literal                tests/b.fastq:fq2 ACAACCTACCGTGACAAAGAAAGTTGTC
check literal-rc             GACAACTTTCTTTTTGTCACGGTAGGTTGT \
                             rc:ACAACCTACCGTGACAAAGAAAGTTGTC
check print-progress         tests/b.fastq tests/b.fastq --print-progress
check copyable-alignment     tests/b.fastq tests/b.fastq \
                             --just-print-copyable-alignment
