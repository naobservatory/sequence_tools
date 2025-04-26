#!/bin/bash
# Compresses an input FASTQ chunk using zstd, streams it to S3, and removes the input chunk on success.
# Usage: ./compress_upload.sh <input_chunk.fastq> <s3_output_path.fastq.zst> [zstd_threads] [zstd_level]

set -euo pipefail # Exit on unhandled non-zero status, undefined variables, and honor failures in pipelines

INPUT_CHUNK="$1"
S3_OUTPUT_PATH="$2"
ZSTD_THREADS="${3:-3}" # for 1M read pair chunks, zstd parallelizes well to 3 threads
ZSTD_LEVEL="${4:-15}" # default high compression level

if [ -z "${INPUT_CHUNK}" ] || [ -z "${S3_OUTPUT_PATH}" ]; then
  echo "Usage: $0 <input_chunk.fastq> <s3_output_path.fastq.zst> [zstd_threads] [zstd_level]" >&2
  exit 1
fi

echo "Starting job for ${INPUT_CHUNK}" >&2
echo "Compressing ${INPUT_CHUNK} to ${S3_OUTPUT_PATH} using ${ZSTD_THREADS} threads and compression level ${ZSTD_LEVEL}" >&2

# Compress, stream to S3 via pipe
zstd -"${ZSTD_LEVEL}" -T"${ZSTD_THREADS}" -c "${INPUT_CHUNK}" | aws s3 cp - "${S3_OUTPUT_PATH}"

# With set -o pipefail enabled, $? will reflect failure of any command in the pipeline,
# not just the aws command
if [ $? -ne 0 ]; then
  echo "ERROR: Compression or S3 upload failed for ${S3_OUTPUT_PATH}" >&2
  # File is not removed on failure to allow for debugging
  exit 1
fi

echo "Upload successful, removing ${INPUT_CHUNK}" >&2
rm "${INPUT_CHUNK}"

echo "Finished job for ${INPUT_CHUNK}" >&2
exit 0
