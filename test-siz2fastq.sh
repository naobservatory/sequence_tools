#!/bin/bash
INVOCATION=""
PREFIX="siz2fastq"

source test-util.sh

check pair1 \
   'cat tests/interleaved.fastq | zstd | ./siz2fastq >(gunzip) /dev/null'
check pair2 \
   'cat tests/interleaved.fastq | zstd | ./siz2fastq /dev/null >(gunzip)'
