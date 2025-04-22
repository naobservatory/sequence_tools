# sequence tools

## align.py

```
usage: align [-h] [--columns N] [--max-dist N] [--min-score N]
             [--print-progress] [--start-pos [1:2]:N] [--moving-average-tsv]
             [--moving-average-width N] [--chart]
             [--just-print-copyable-alignment]
             [rc:]ACTG|path[:id] [:rc]ACTG|path[:id]

Align sequences and show them vertically interleaved

positional arguments:
  [rc:]ACTG|path[:id]   First input. Either a literal sequence or a path to a
                        fasta or fastq file. If a path, optionally include a
                        colon-separated id to refer to identify a specific
                        record in the file.
  [:rc]ACTG|path[:id]   Same as first input.

options:
  -h, --help            show this help message and exit
  --columns N           How many columns to wrap at. If unspecified,
                        autodetects.
  --max-dist N          How hard to try to avoid over-alignment when blocks
                        have been substituted out.
  --min-score N         Minimum score of alignment to print.
  --print-progress      Mark lines with how far along the genomes they
                        represent.
  --start-pos [1:2]:N   Start position, as bases into either sequence 1 or 1
  --moving-average-tsv  Instead of printing bases, output a tsv of how well
                        these genomes align.
  --moving-average-width N
                        Standard deviation of the moving average. Ignored
                        unless --moving-average-tsv.
  --chart               visualize the alignment with some squiggles
  --just-print-copyable-alignment
                        Ignore any other settings and just use the default
                        BioPython alignment pretty printing
```

## siz2fastq

```
Usage: ./siz2fastq [-z] out1.fastq[.gz] out2.fastq[.gz]
Options:
  -z    Compress output with gzip
```

Convert zstd-compressed interleaved fastq to paired (optionally
gzip-compressed) fastq.

Unlike full fastq, hard wrapping is not permitted: this simply takes the first
four lines and puts them in `output_1.fastq.gz`, then the next four and puts
them in `output_2.fastq.gz`, and then keeps alternating.

Compatible with streaming. For example:

```
aws s3 cp s3://bucket/foo.fastq.zstd | \
    siz2fastq -z >(aws s3 cp - s3://bucket/foo_1.fastq.gz) \
                 >(aws s3 cp - s3://bucket/foo_2.fastq.gz) \
```

## SIZer

Convert paired short-read FASTQ files to split, interleaved, zstd-compressed (SIZ) format.
Typically used within a data processing pipeline and with compress-upload.sh.

Presently, input FASTQ files with hard wrapping are not supported.
```
Usage: ./SIZer [-p <local_prefix>] <reads_per_chunk> <r1_fastq> <r2_fastq> <s3_output_prefix> <compress_script_path> <max_concurrent_jobs> <zstd_threads>
Options:
  -p <local_prefix>  Prefix for temporary chunk files (default: work_chunk)
Arguments:
  reads_per_chunk      Max read pairs per output chunk
  r1_fastq             Path to R1 FASTQ file (can be .gz)
  r2_fastq             Path to R2 FASTQ file (can be .gz)
  s3_output_prefix     S3 prefix for output files
  compress_script_path Path to compression/upload script
  max_concurrent_jobs  Max simultaneous background jobs
  zstd_threads         Number of zstd threads for script

Example:
  ./SIZer -p work/sampleA 1000000 r1.fq.gz r2.fq.gz s3://mybucket/siz/sampleA ./compress_upload.sh 16 3
```
