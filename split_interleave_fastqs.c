#include <dirent.h>
#include <errno.h>
#include <getopt.h>
#include <libgen.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define MAX_SUFFIX_LEN 24  // _chunkNNNNNN.fastq
#define MAX_FILENAME_LEN 4096
#define DEFAULT_READ_PAIRS_PER_CHUNK 1000000

#ifdef ENABLE_DEBUG
  #define DEBUG(fmt, ...) fprintf(stderr, "Debug: " fmt "\n", ##__VA_ARGS__)
#else
  #define DEBUG(fmt, ...) ((void)0)
#endif

// Buffers for a 4-line FASTQ record
typedef struct {
  char *lines[4];
  size_t sizes[4];
} FastqRecordBuffers;

void print_usage_and_exit(const char *program_name);
FILE *open_input_or_die(const char *path);
void close_file(FILE *f, const char *path);
int read_fastq_record(FILE *f, FastqRecordBuffers *buffers, int input_idx);
void cleanup_buffers(FastqRecordBuffers *buffers);
void perror_and_exit(const char *context, const char *details);
int count_chunk_files(const char *dir_path, const char *base_name);
void wait_for_file_slot(const char *dir_path, const char *base_name, int max_files);

int main(int argc, char **argv) {
  char *output_prefix = NULL;
  long max_files = -1;
  long read_pairs_per_chunk = DEFAULT_READ_PAIRS_PER_CHUNK;
  int opt;

  // Parse command line arguments
  while ((opt = getopt(argc, argv, "p:n:c:")) != -1) {
    switch (opt) {
      case 'p':
        output_prefix = optarg;
        break;
      case 'n':
        errno = 0;
        max_files = strtol(optarg, NULL, 10);
        if (max_files <= 0) {
          fprintf(stderr, "Error: Invalid -n value '%s'. Must be positive integer.\n",
                  optarg);
          print_usage_and_exit(argv[0]);
        }
        break;
      case 'c':
        errno = 0;
        read_pairs_per_chunk = strtol(optarg, NULL, 10);
        if (read_pairs_per_chunk <= 0) {
          fprintf(stderr, "Error: Invalid -c value '%s'. Must be positive integer.\n",
                  optarg);
          print_usage_and_exit(argv[0]);
        }
        break;
      default:
        print_usage_and_exit(argv[0]);
    }
  }

  // Validate arguments
  if (output_prefix == NULL || max_files == -1) {
    fprintf(stderr, "Error: Missing required argument(s) -p or -n.\n");
    print_usage_and_exit(argv[0]);
  }
  if (optind != argc - 2) {
    fprintf(stderr,
            "Error: Expected exactly two positional arguments: <r1_fastq> <r2_fastq>\n");
    print_usage_and_exit(argv[0]);
  }

  char *r1_path = argv[optind];
  char *r2_path = argv[optind + 1];

  // Extract directory and basename from output_prefix
  // Need separate copies because dirname/basename can modify their inputs
  char *prefix_copy_dir = strdup(output_prefix);
  char *prefix_copy_base = strdup(output_prefix);
  if (!prefix_copy_dir || !prefix_copy_base) {
    free(prefix_copy_dir);
    free(prefix_copy_base);
    perror_and_exit("strdup prefix copy", output_prefix);
  }

  char *dir_path = dirname(prefix_copy_dir);
  char *base_name = basename(prefix_copy_base);

  // Initialize resources
  FILE *f1 = open_input_or_die(r1_path);
  FILE *f2 = open_input_or_die(r2_path);
  FILE *out_chunk_fp = NULL;
  FastqRecordBuffers r1_buffers = {{NULL}, {0}};
  FastqRecordBuffers r2_buffers = {{NULL}, {0}};
  char chunk_suffix[MAX_SUFFIX_LEN];
  char chunk_filename[MAX_FILENAME_LEN];
  int chunk_index = 0;
  long current_chunk_reads = 0;
  long total_read_pairs = 0;
  int exit_code = 0;

  DEBUG("Starting: prefix='%s', max_files=%ld, reads/chunk=%ld\n", output_prefix,
        max_files, read_pairs_per_chunk);

  // Main processing loop - read and interleave FASTQ records
  while (true) {
    bool r1_eof = (read_fastq_record(f1, &r1_buffers, 1) == 0);
    bool r2_eof = (read_fastq_record(f2, &r2_buffers, 2) == 0);

    if (r1_eof != r2_eof) {
      fprintf(stderr, "Error: Unequal reads in R1 ('%s') and R2 ('%s').\n", r1_path,
              r2_path);
      exit_code = 1;
      break;
    }

    if (r1_eof) break;

    if (out_chunk_fp == NULL) {
      // Wait for a slot and open new chunk file
      wait_for_file_slot(dir_path, base_name, max_files);
      snprintf(chunk_suffix, MAX_SUFFIX_LEN, "_chunk%06d.fastq", chunk_index);
      snprintf(chunk_filename, MAX_FILENAME_LEN, "%s/%s%s", dir_path,
               base_name, chunk_suffix);
      out_chunk_fp = fopen(chunk_filename, "w");
      if (!out_chunk_fp) perror_and_exit("fopen output chunk", chunk_filename);
      DEBUG("Opened chunk '%s'\n", chunk_filename);
      current_chunk_reads = 0;
    }

    // Write interleaved records
    for (int i = 0; i < 4; ++i) {
      if (fputs(r1_buffers.lines[i], out_chunk_fp) == EOF)
        perror_and_exit("fputs to chunk", chunk_filename);
    }
    for (int i = 0; i < 4; ++i) {
      if (fputs(r2_buffers.lines[i], out_chunk_fp) == EOF)
        perror_and_exit("fputs to chunk", chunk_filename);
    }

    total_read_pairs++;
    current_chunk_reads++;

    // Close chunk when full
    if (current_chunk_reads >= read_pairs_per_chunk) {
      if (fclose(out_chunk_fp) != 0) perror_and_exit("fclose chunk", chunk_filename);
      DEBUG("Closed chunk '%s' (%ld pairs)\n", chunk_filename, current_chunk_reads);
      printf("%s\n", chunk_suffix);  // Output chunk suffix to stdout, not stderr
      fflush(stdout);
      out_chunk_fp = NULL;
      chunk_index++;
    }
  }

  // Handle final chunk
  if (exit_code == 0 && out_chunk_fp != NULL) {
    if (fclose(out_chunk_fp) != 0) perror_and_exit("fclose final chunk", chunk_filename);
    DEBUG("Closed final chunk '%s' (%ld pairs)\n", chunk_filename, current_chunk_reads);
    printf("%s\n", chunk_suffix);
    fflush(stdout);
    out_chunk_fp = NULL;
  }

  // Cleanup
  fprintf(stderr, "Finished. Processed %ld total pairs.\n", total_read_pairs);
  close_file(f1, r1_path);
  close_file(f2, r2_path);
  cleanup_buffers(&r1_buffers);
  cleanup_buffers(&r2_buffers);
  if (out_chunk_fp != NULL) {
    DEBUG("Closing incomplete chunk '%s' due to error.\n", chunk_filename);
    fclose(out_chunk_fp);
  }

  free(prefix_copy_dir);
  free(prefix_copy_base);

  DEBUG("Exiting with code %d.\n", exit_code);
  return exit_code;
}

void print_usage_and_exit(const char *program_name) {
  fprintf(stderr,
          "Usage: %s -p <prefix> -n <max_files> [-c <reads>] <r1_fastq> <r2_fastq>\n",
          program_name);
  fprintf(stderr, "\nOptions:\n");
  fprintf(stderr, "  -p <prefix>    (Required) Path prefix for output chunk files; ");
  fprintf(stderr, "a memory-backed location is recommended.\n");
  fprintf(stderr, "                 (e.g., /dev/shm/myfastq -> ");
  fprintf(stderr, "/dev/shm/myfastq_chunkNNNNNN.fastq)\n");
  fprintf(stderr,
          "  -n <max_files> (Required) Max number of chunk files allowed concurrently.\n");
  fprintf(stderr, "  -c <reads>     (Optional) Target read pairs per chunk (default: %d).\n",
          DEFAULT_READ_PAIRS_PER_CHUNK);
  fprintf(stderr, "\nArguments:\n");
  fprintf(stderr, "  r1_fastq      Path to R1 FASTQ file.\n");
  fprintf(stderr, "  r2_fastq      Path to R2 FASTQ file.\n");
  fprintf(stderr, "\nExample:\n");
  fprintf(stderr, "  %s -p /tmp/sampleA -n 16 -c 1000000 sA_R1.fq sA_R2.fq\n",
          program_name);
  exit(1);
}

void perror_and_exit(const char *context, const char *details) {
  fprintf(stderr, "Error: %s (%s): ", context, details);
  perror(NULL);
  exit(1);
}

FILE *open_input_or_die(const char *path) {
  FILE *f = fopen(path, "r");
  if (!f) perror_and_exit("fopen input", path);
  DEBUG("Opened input '%s'.\n", path);
  return f;
}

void close_file(FILE *f, const char *path) {
  if (f == NULL) return;
  if (fclose(f) != 0) {
    fprintf(stderr, "Warning: fclose failed for '%s': %s\n", path, strerror(errno));
  } else {
    DEBUG("Closed '%s'.\n", path);
  }
}

// Read 4 lines into buffers with no parsing: we just treat every 4 lines as a record.
// input_idx is 1 or 2, to indicate which input was problematic in the case of an error.
int read_fastq_record(FILE *f, FastqRecordBuffers *buffers, int input_idx) {
  for (int i = 0; i < 4; ++i) {
    errno = 0;
    ssize_t len = getline(&(buffers->lines[i]), &(buffers->sizes[i]), f);
    if (len == -1) {
      if (i == 0 && feof(f)) return 0;  // Clean EOF before record
      if (feof(f))
        fprintf(stderr, "Error: Unexpected EOF in input %d.\n", input_idx);
      else {
        fprintf(stderr, "Error reading line from input %d: ", input_idx);
        perror(NULL);
      }
      exit(1);
    }
  }
  return 1;
}

void cleanup_buffers(FastqRecordBuffers *buffers) {
  for (int i = 0; i < 4; ++i) {
    free(buffers->lines[i]);
    buffers->lines[i] = NULL;
  }
}

// Returns count of files in dir_path matching "<base_name>_chunk*.fastq" pattern
int count_chunk_files(const char *dir_path, const char *base_name) {
  DIR *dir = opendir(dir_path);
  if (!dir) {
    if (errno == ENOENT) return 0;  // Dir doesn't exist = 0 files
    perror_and_exit("opendir in count_chunk_files", dir_path);
  }

  struct dirent *entry;
  int count = 0;
  char prefix[MAX_FILENAME_LEN];
  snprintf(prefix, MAX_FILENAME_LEN, "%s_chunk", base_name);
  size_t prefix_len = strlen(prefix);
  const char *fastq_ext = ".fastq";
  size_t fastq_ext_len = strlen(fastq_ext);

  while ((entry = readdir(dir)) != NULL) {
    const char *name = entry->d_name;
    size_t name_len = strlen(name);

    // Count files matching pattern: <base_name>_chunk*.fastq
    if (name_len >= prefix_len + fastq_ext_len && strncmp(name, prefix, prefix_len) == 0 &&
        strcmp(name + name_len - fastq_ext_len, fastq_ext) == 0) {
      count++;
    }
  }

  if (closedir(dir) == -1) {
    fprintf(stderr, "Warning: closedir failed for '%s': %s\n", dir_path, strerror(errno));
  }

  return count;
}

// Blocks until number of chunk files is less than max_files
void wait_for_file_slot(const char *dir_path, const char *base_name, int max_files) {
  bool reported_waiting = false;
  while (true) {
    int current_files = count_chunk_files(dir_path, base_name);

    if (current_files < max_files) {
      if (reported_waiting) {
        DEBUG("Slot available (%d/%d files). Resuming.\n", current_files,
                max_files);
      }
      return;
    }

    if (!reported_waiting) {
      DEBUG("Max files limit (%d/%d) reached. Waiting...\n", current_files,
              max_files);
      reported_waiting = true;
    }

    usleep(50000);  // Wait 0.05 seconds
  }
}

