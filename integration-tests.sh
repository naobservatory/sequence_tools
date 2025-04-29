#!/bin/bash
set -e

EXPECTED_HASH="a000ca9241a286441bc46bfdb5c12723" 

echo "Running sizer command..."
./sizer.sh <(aws s3 cp s3://nao-testing/sequence-tools/sizer/raw/test01_1.fastq.gz - | gunzip) \
           <(aws s3 cp s3://nao-testing/sequence-tools/sizer/raw/test01_2.fastq.gz - | gunzip) \
           sizer_integration_test \
           s3://nao-testing/sequence-tools/sizer/siz/test01

# Wait briefly for S3 consistency
sleep 1

# md5sum will give us a trailing " -" since we provide stdin input, so we cut that off
if diff <(aws s3 cp s3://nao-testing/read-sizer/siz/test01_chunk000000.fastq.zst - | md5sum | cut -d' ' -f1) \
        <(echo "$EXPECTED_HASH") >/dev/null; then
    echo "INTEGRATION TEST PASSED: Hashes match"
    exit 0
else
    echo "INTEGRATION TEST FAILED: Hashes do not match"
    exit 1
fi
