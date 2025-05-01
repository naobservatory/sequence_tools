#!/bin/bash
# Wrapper script for processing FASTQ pairs to Split-Interleaved-Zstd (SIZ) files in S3
# Usage: ./fastq_to_siz.sh [options] <r1_fastq> <r2_fastq> <local_prefix> <s3_output_prefix>

set -euo pipefail

# Default values
TEMPDIR="/tmp/split_interleave"
PARALLEL_TASKS=$(( $(nproc) / 3 ))  # Floor of CPU count / 3
ZSTD_THREADS=3
ZSTD_LEVEL=15
READ_PAIRS_PER_CHUNK=1000000
SPLIT_INTERLEAVE_BIN="./split_interleave_fastqs"
COMPRESS_UPLOAD_BIN="./compress_upload.sh"

print_usage_and_exit() {
    echo "Usage: $0 [options] <r1_fastq> <r2_fastq> <prefix> <s3_output_prefix>"
    echo ""
    echo "Arguments:"
    echo "  <r1_fastq>          Forward read FASTQ file or stream"
    echo "  <r2_fastq>          Reverse read FASTQ file or stream"
    echo "  <prefix>            Filename prefix for temporary files"
    echo "  <s3_output_prefix>  S3 path prefix for output SIZ files (e.g., s3://bucket/path/sample)"
    echo ""
    echo "Options:"
    echo "  -d, --tempdir DIR         Local directory for temporary split-interleaved files. Should be memory backed. (default: $TEMPDIR)"
    echo "  -j, --jobs JOBS           Number of parallel compress-upload tasks (default: $PARALLEL_TASKS)"
    echo "  -t, --threads THREADS     ZSTD threads per compress-upload task (default: $ZSTD_THREADS)"
    echo "  -l, --level LEVEL         ZSTD compression level (default: $ZSTD_LEVEL)"
    echo "  -c, --chunk-size SIZE     Read pairs per SIZ file (default: $READ_PAIRS_PER_CHUNK)"
    echo "  -s, --split-bin PATH      Path to split_interleave_fastqs binary (default: $SPLIT_INTERLEAVE_BIN)"
    echo "  -u, --upload-bin PATH     Path to compress_upload.sh script (default: $COMPRESS_UPLOAD_BIN)"
    echo "  -h, --help                Display this help message"
    echo ""
    exit 1
}

# Parse arguments
POSITIONAL_ARGS=()
while [[ $# -gt 0 ]]; do
  case $1 in
    -d|--tempdir)
      TEMPDIR="$2"
      shift 2
      ;;
    -j|--jobs)
      PARALLEL_TASKS="$2"
      shift 2
      ;;
    -t|--threads)
      ZSTD_THREADS="$2"
      shift 2
      ;;
    -l|--level)
      ZSTD_LEVEL="$2"
      shift 2
      ;;
    -c|--chunk-size)
      READ_PAIRS_PER_CHUNK="$2"
      shift 2
      ;;
    -s|--split-bin)
      SPLIT_INTERLEAVE_BIN="$2"
      shift 2
      ;;
    -u|--upload-bin)
      COMPRESS_UPLOAD_BIN="$2"
      shift 2
      ;;
    -h|--help)
      print_usage_and_exit
      ;;
    -*|--*)
      echo "Error: Unknown option $1"
      print_usage_and_exit
      ;;
    *)
      POSITIONAL_ARGS+=("$1")
      shift
      ;;
  esac
done

# Restore positional parameters
set -- "${POSITIONAL_ARGS[@]}"

if [ $# -ne 4 ]; then
    echo "Error: Expected exactly 4 positional arguments."
    print_usage_and_exit
fi

R1_FASTQ="$1"
R2_FASTQ="$2"
PREFIX="$3"
S3_OUTPUT_PREFIX="$4"

if [ ! -x "$SPLIT_INTERLEAVE_BIN" ]; then
    echo "Error: Split-interleave binary $SPLIT_INTERLEAVE_BIN does not exist or is not executable."
    exit 1
fi
if [ ! -x "$COMPRESS_UPLOAD_BIN" ]; then
    echo "Error: Compress-upload script $COMPRESS_UPLOAD_BIN does not exist or is not executable."
    exit 1
fi

# Check for existing files with the PREFIX
if ls "${TEMPDIR}/${PREFIX}"_chunk*.fastq 1> /dev/null 2>&1; then
    echo "Error: Files matching ${TEMPDIR}/${PREFIX}_chunk*.fastq already exist. Please remove them or use a different prefix."
    exit 1
fi

echo "Starting SIZ pipeline with the following parameters:"
echo "  Forward reads:       $R1_FASTQ"
echo "  Reverse reads:       $R2_FASTQ"
echo "  Local prefix:        $PREFIX"
echo "  S3 output prefix:    $S3_OUTPUT_PREFIX"
echo "  Parallel tasks:      $PARALLEL_TASKS"
echo "  ZSTD threads:        $ZSTD_THREADS"
echo "  ZSTD level:          $ZSTD_LEVEL"
echo "  Read pairs/chunk:    $READ_PAIRS_PER_CHUNK"
echo ""

# Run the pipeline, using xargs for compress-upload parallelization
mkdir -p "$TEMPDIR"
MAX_FILES=$((PARALLEL_TASKS + 2)) # keep each compression task busy
"$SPLIT_INTERLEAVE_BIN" -d "$TEMPDIR" -p "$PREFIX" -n "$MAX_FILES" -c "$READ_PAIRS_PER_CHUNK" "$R1_FASTQ" "$R2_FASTQ" | \
    xargs -P "$PARALLEL_TASKS" -I {} "$COMPRESS_UPLOAD_BIN" "${TEMPDIR}/${PREFIX}{}" "${S3_OUTPUT_PREFIX}{}.zst" "$ZSTD_THREADS" "$ZSTD_LEVEL"

echo "Pipeline completed successfully."
