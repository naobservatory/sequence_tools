#!/bin/bash
INVOCATION=""
PREFIX="split_interleave_fastqs"

source test-util.sh

# Test 1: Interleave b.fastq with itself, making 2 chunks with 2 read pairs per chunk
check interleave-chunk \
   'mkdir -p /tmp/test-split-interleave && rm -f /tmp/test-split-interleave/test_chunk* && ./split_interleave_fastqs -p /tmp/test-split-interleave/test -n 2 -c 2 tests/b.fastq tests/b.fastq 2>/dev/null; rm -rf /tmp/test-split-interleave'

# Test 2: Interleave with -c 3 -n 1 and cat the output to see the actual interleaved content
check interleave-cat-content \
   'mkdir -p /tmp/test-split-interleave && rm -f /tmp/test-split-interleave/test_chunk* && ./split_interleave_fastqs -p /tmp/test-split-interleave/test -n 1 -c 3 tests/b.fastq tests/b.fastq 2>/dev/null > /tmp/test-split-interleave/file_list && cat $(cat /tmp/test-split-interleave/file_list) && rm -rf /tmp/test-split-interleave'
