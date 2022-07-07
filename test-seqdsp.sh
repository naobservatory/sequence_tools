#!/bin/bash
INVOCATION="./seqdsp.py"
PREFIX="seqdsp"
source test-util.sh

check display-fasta          tests/a.fasta
check fasta-by-id            tests/a.fasta:seq123
check seq-matches            tests/a.fasta --seq-matches GGGCCTATATAGGATC
check seq-matches-rc         tests/a.fasta --seq-matches rcGATGCTATTTAGGCCC
check seq-matches-color      tests/a.fasta --seq-matches GGGCCTATATAGGATC:blue
check seq-matches-rc-color   tests/a.fasta --seq-matches rcGATGCTATTTAGGCCC:blue
check seq-matches-multiple   tests/a.fasta --seq-matches GGGCCTATAT \
                                           --seq-matches AAAATCGC:green
check seq-matches-different  tests/a.fasta --seq-matches GGGCCTATATAGGATC:blue \
                                           --seq-matches TTTTATGATTTATACTATACG
check display-fastq          tests/b.fastq
check show-quality           tests/b.fastq --show-quality
check colorize-quality       tests/b.fastq --show-quality --colorize-quality
check max-quality            tests/b.fastq --show-quality --colorize-quality \
                                           --max-quality '$'
check showquality-when-highlighted \
                             tests/b.fastq --show-quality-when-highlighted \
                                           --seq-matches CCGATAGCTCTACCGATTGAAA
check skip-lines             tests/b.fastq --show-quality-when-highlighted \
                                           --skip-lines \
                                           --seq-matches CCGATAGCTCTACCGATTGAAA
                             
