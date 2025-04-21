#!/bin/bash
INPUT_CHUNK="$1"
TARGET_S3_PATH="$2"
ZSTD_THREADS="$3"

echo "Script Args: Input=${INPUT_CHUNK} Output=${TARGET_S3_PATH} Threads=${ZSTD_THREADS}"
echo "--- Chunk Content Start ---"
cat "${INPUT_CHUNK}"
echo "--- Chunk Content End ---"

# Clean up the uncompressed input chunk
rm "${INPUT_CHUNK}"

exit 0

