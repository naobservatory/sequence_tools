#!/bin/bash

INVOCATION="./SIZer"
PREFIX="SIZer/" # for naming gold files

source test-util.sh

# check <test_name> <reads/chunk> <r1> <r2> <s3_prefix> <script> <jobs> <threads>
check combined_chunking \
    2 \
    tests/SIZer/input_1.fastq \
    tests/SIZer/input_2.fastq \
    s3://dummy/s3_out \
    tests/SIZer/dummy_compress_upload.sh \
    1 \
    1

echo "SIZer tests finished."
