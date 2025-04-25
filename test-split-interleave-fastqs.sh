#!/bin/bash
INVOCATION=""
PREFIX="split_interleave_fastqs"

source test-util.sh

LOCAL_DIR="/tmp/test-split-interleave"
LOCAL_PREFIX="test_sample"
LOCAL_BASE_PATH="${LOCAL_DIR}/${LOCAL_PREFIX}"

# Interleave b.fastq with itself, making 2 chunks with 2 read pairs per chunk
check interleave-chunk \
   "mkdir -p ${LOCAL_DIR} && rm -f ${LOCAL_BASE_PATH}_chunk* && ./split_interleave_fastqs -d ${LOCAL_DIR} -n 2 -c 2 tests/b.fastq tests/b.fastq 2>/dev/null; rm -rf ${LOCAL_DIR}"

# Same, but with optional prefix argument
check interleave-chunk-with-p \
   "mkdir -p ${LOCAL_DIR} && rm -f ${LOCAL_BASE_PATH}_chunk* && ./split_interleave_fastqs -d ${LOCAL_DIR} -n 2 -c 2 -p ${LOCAL_PREFIX} tests/b.fastq tests/b.fastq 2>/dev/null; rm -rf ${LOCAL_DIR}"

# Test 2: Interleave with -c 3 -n 1 and cat the output to see the actual interleaved content
check interleave-cat-content \
   "mkdir -p ${LOCAL_DIR} && rm -f ${LOCAL_BASE_PATH}_chunk* && ./split_interleave_fastqs -d ${LOCAL_DIR} -n 1 -c 3 -p ${LOCAL_PREFIX} tests/b.fastq tests/b.fastq 2>/dev/null | xargs -I {} cat ${LOCAL_BASE_PATH}{} && rm -rf ${LOCAL_DIR}"
