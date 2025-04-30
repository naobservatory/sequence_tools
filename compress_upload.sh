#!/bin/bash
# Compresses an input FASTQ chunk using zstd, streams it to S3 or a local file, and
# removes the input chunk on success.

set -euo pipefail # Exit on unhandled non-zero status, undefined variables, and honor failures in pipelines

INPUT_CHUNK="$1"
OUTPUT_PATH="$2"
ZSTD_THREADS="${3:-3}" # for 1M read pair chunks, zstd parallelizes well to 3 threads
ZSTD_LEVEL="${4:-15}" # default high compression level

if [ -z "${INPUT_CHUNK}" ] || [ -z "${OUTPUT_PATH}" ]; then
  echo "Usage: $0 <input_chunk.fastq> <output_path.fastq.zst> [zstd_threads] [zstd_level]" >&2
  exit 1
fi

echo "Starting job for ${INPUT_CHUNK}" >&2
echo "Compressing ${INPUT_CHUNK} to ${OUTPUT_PATH} using ${ZSTD_THREADS} threads and compression level ${ZSTD_LEVEL}" >&2

# Compress and stream; copy mechanism determined by destination type
if [[ "${OUTPUT_PATH}" == s3://* ]]; then
  zstd -"${ZSTD_LEVEL}" -T"${ZSTD_THREADS}" -c "${INPUT_CHUNK}" | aws s3 cp - "${OUTPUT_PATH}"
else
  # Local path (including /dev/null)
  zstd -"${ZSTD_LEVEL}" -T"${ZSTD_THREADS}" -c "${INPUT_CHUNK}" > "${OUTPUT_PATH}"
fi

echo "Upload successful, removing ${INPUT_CHUNK}" >&2
rm "${INPUT_CHUNK}"

echo "Finished job for ${INPUT_CHUNK}" >&2
exit 0
