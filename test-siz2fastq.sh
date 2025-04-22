#!/bin/bash
INVOCATION=""
PREFIX="siz2fastq"

source test-util.sh

check pair1 \
   'cat tests/interleaved.fastq | zstd | ./siz2fastq >(cat) /dev/null'
check pair2 \
   'cat tests/interleaved.fastq | zstd | ./siz2fastq /dev/null >(cat)'
check pair1-gz \
   'cat tests/interleaved.fastq | zstd | ./siz2fastq -z >(gunzip) /dev/null'
check pair2-gz \
   'cat tests/interleaved.fastq | zstd | ./siz2fastq -z /dev/null >(gunzip)'
