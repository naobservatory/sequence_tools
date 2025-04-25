#!/bin/bash
INVOCATION=""
PREFIX="split_interleave_fastqs"

source test-util.sh

LOCAL_DIR="/tmp/test-split-interleave"
LOCAL_PREFIX="${LOCAL_DIR}/sample"

# Test 1: Interleave b.fastq with itself, making 2 chunks with 2 read pairs per chunk
check interleave-chunk \
   "mkdir -p ${LOCAL_DIR} && rm -f ${LOCAL_PREFIX}_chunk* && ./split_interleave_fastqs -p ${LOCAL_PREFIX} -n 2 -c 2 tests/b.fastq tests/b.fastq 2>/dev/null; rm -rf ${LOCAL_DIR}"

# Test 2: Interleave with -c 3 -n 1 and cat the output to see the actual interleaved content
check interleave-cat-content \
   "mkdir -p ${LOCAL_DIR} && rm -f ${LOCAL_PREFIX}_chunk* && ./split_interleave_fastqs -p ${LOCAL_PREFIX} -n 1 -c 3 tests/b.fastq tests/b.fastq 2>/dev/null | xargs -I {} cat ${LOCAL_PREFIX}{} && rm -rf ${LOCAL_DIR}"
