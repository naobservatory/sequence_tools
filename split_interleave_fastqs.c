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
#define DEFAULT_TEMP_PREFIX "sample"
#define IO_BUF_SIZE 8192

// State for one input stream (R1 or R2)
typedef struct {
  FILE *f;
  const char *path;
  char *buf;
  size_t buf_size;
  size_t pos;
  size_t valid;
  bool eof_flag;  // indicates clean EOF was reached by fread
} InputState;

// State for the shared output buffer
typedef struct {
  FILE *f;  // NULL if no chunk open
  char *buf;
  size_t buf_size;
  size_t pos;
} OutputState;

#ifdef ENABLE_DEBUG
#define DEBUG(fmt, ...) fprintf(stderr, "Debug: " fmt, ##__VA_ARGS__)
#else
#define DEBUG(fmt, ...) ((void)0)
#endif

void print_usage_and_exit(const char *program_name);
FILE *open_input_or_die(const char *path);
void close_file(FILE *f, const char *path);
void perror_and_exit(const char *context, const char *details);
int count_chunk_files(const char *dir_path, const char *base_name);
void wait_for_file_slot(const char *dir_path, const char *base_name, int max_files);
int copy_record(InputState *input, OutputState *output, const char *out_path);

int main(int argc, char **argv) {
  char *output_dir = NULL;
  char *output_prefix = DEFAULT_TEMP_PREFIX;
  long max_files = -1;
  long read_pairs_per_chunk = DEFAULT_READ_PAIRS_PER_CHUNK;
  int opt;

  // Parse command line arguments
  while ((opt = getopt(argc, argv, "d:n:c:p:")) != -1) {
    switch (opt) {
      case 'd':
        output_dir = optarg;
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
      case 'p':
        output_prefix = optarg;
        break;
      default:
        print_usage_and_exit(argv[0]);
    }
  }

  // Validate arguments
  if (output_dir == NULL || max_files == -1) {
    fprintf(stderr, "Error: Missing required argument(s) -d or -n.\n");
    print_usage_and_exit(argv[0]);
  }
  if (optind != argc - 2) {
    fprintf(stderr,
            "Error: Expected exactly two positional arguments: <r1_fastq> <r2_fastq>\n");
    print_usage_and_exit(argv[0]);
  }

  char *r1_path = argv[optind];
  char *r2_path = argv[optind + 1];

  // Initialize resources
  // inputs
  char r1_in_buf[IO_BUF_SIZE];
  char r2_in_buf[IO_BUF_SIZE];

  InputState r1_state = {
      open_input_or_die(r1_path), r1_path, r1_in_buf, IO_BUF_SIZE, 0, 0, false};
  InputState r2_state = {
      open_input_or_die(r2_path), r2_path, r2_in_buf, IO_BUF_SIZE, 0, 0, false};
  // output
  char out_buf[IO_BUF_SIZE];
  OutputState out_state = {NULL, out_buf, IO_BUF_SIZE, 0};
  // state
  char chunk_suffix[MAX_SUFFIX_LEN];
  char chunk_filename[MAX_FILENAME_LEN];
  int chunk_index = 0;
  long current_chunk_reads = 0;
  long total_read_pairs = 0;
  int exit_code = 0;

  DEBUG("Starting: prefix='%s', max_files=%ld, reads/chunk=%ld\n", output_prefix, max_files,
        read_pairs_per_chunk);

  // Main processing loop - read and interleave FASTQ records
  while (true) {
    if (out_state.f == NULL) {
      // Wait for a slot and open new chunk file
      wait_for_file_slot(output_dir, output_prefix, max_files);
      snprintf(chunk_suffix, MAX_SUFFIX_LEN, "_chunk%06d.fastq", chunk_index);
      snprintf(chunk_filename, MAX_FILENAME_LEN, "%s/%s%s", output_dir, output_prefix,
               chunk_suffix);
      out_state.f = fopen(chunk_filename, "w");
      if (!out_state.f) perror_and_exit("fopen output chunk", chunk_filename);
      DEBUG("Opened chunk '%s'\n", chunk_filename);
      current_chunk_reads = 0;
    }

    // copy a fastq record from each input
    int status_r2;
    int status_r1 = copy_record(&r1_state, &out_state, chunk_filename);
    if (status_r1 >= 0) {  // Don't process R2 if R1 had error/unexpected EOF
      status_r2 = copy_record(&r2_state, &out_state, chunk_filename);
    } else {
      status_r2 = -1;  // treat R2 as failed if R1 failed: we didn't interleave a full pair
    }

    // are we done copying?
    if (status_r1 == 0 && status_r2 == 0) {  // Both ended cleanly at record boundary
      DEBUG("Clean EOF detected on both inputs.\n");
      break;  // normal loop exit
    }
    if (status_r1 == -1 || status_r2 == -1) {  // Error or unexpected EOF in either
      fprintf(stderr, "Exiting due to error or unexpected EOF.\n");
      exit_code = 1;
      break;
    }
    if (status_r1 != status_r2) {
      fprintf(stderr,
              "Error: Unequal reads detected between '%s' and '%s'. One file ended before "
              "the other at record boundary.\n",
              r1_path, r2_path);
      exit_code = 1;
      break;
    }

    total_read_pairs++;
    current_chunk_reads++;

    // Close chunk when full
    if (current_chunk_reads >= read_pairs_per_chunk) {
      // flush and reset output buffer
      if (out_state.pos > 0) {
        // size of 1 for buffer of characters
        if (fwrite(out_state.buf, 1, out_state.pos, out_state.f) != out_state.pos) {
          perror_and_exit("fwrite error for full buffer in copy_record", chunk_filename);
        }
        out_state.pos = 0;
      }
      // close file
      if (fclose(out_state.f) != 0) perror_and_exit("fclose chunk", chunk_filename);
      DEBUG("Closed chunk '%s' (%ld pairs)\n", chunk_filename, current_chunk_reads);
      printf("%s\n", chunk_suffix);  // Output chunk suffix to stdout, not stderr
      fflush(stdout);
      out_state.f = NULL;
      chunk_index++;
    }
  }  // end main while(true) copy loop

  // Handle final chunk
  if (exit_code == 0 && out_state.f != NULL) {
    // flush output buffer
    if (out_state.pos > 0) {
      if (fwrite(out_state.buf, 1, out_state.pos, out_state.f) != out_state.pos) {
           perror_and_exit("fwrite error flushing final output buffer", chunk_filename);
      }
      out_state.pos = 0;
    }
    if (fclose(out_state.f) != 0) perror_and_exit("fclose final chunk", chunk_filename);
    DEBUG("Closed final chunk '%s' (%ld pairs)\n", chunk_filename, current_chunk_reads);
    printf("%s\n", chunk_suffix);
    fflush(stdout);
    out_state.f = NULL;
  }

  // Cleanup
  fprintf(stderr, "Finished. Processed %ld total pairs.\n", total_read_pairs);
  close_file(r1_state.f, r1_path);
  close_file(r2_state.f, r2_path);
  if (out_state.f != NULL) {
    DEBUG("Closing incomplete chunk '%s' due to error.\n", chunk_filename);
    fclose(out_state.f);
  }

  DEBUG("Exiting with code %d.\n", exit_code);
  return exit_code;
}

void print_usage_and_exit(const char *program_name) {
  fprintf(stderr,
          "Usage: %s -d <temp_dir> -n <max_files> [-c <reads>] [-p <prefix>] <r1_fastq> "
          "<r2_fastq>\n",
          program_name);
  fprintf(stderr, "\nOptions:\n");
  fprintf(stderr, "  -d <temp_dir>  (Required) Directory for output chunk files; ");
  fprintf(stderr, "a memory-backed location is recommended.\n");
  fprintf(stderr,
          "  -n <max_files> (Required) Max number of chunk files allowed concurrently.\n");
  fprintf(stderr,
          "  -c <reads>     (Optional) Target read pairs per chunk (default: %d).\n",
          DEFAULT_READ_PAIRS_PER_CHUNK);
  fprintf(stderr,
          "  -p <prefix>    (Optional) Prefix for temp files in <temp_dir> (default: %s)",
          DEFAULT_TEMP_PREFIX);
  fprintf(stderr, "\nArguments:\n");
  fprintf(stderr, "  r1_fastq      Path to R1 FASTQ file.\n");
  fprintf(stderr, "  r2_fastq      Path to R2 FASTQ file.\n");
  fprintf(stderr, "\nExample:\n");
  fprintf(stderr, "  %s -d /tmp -n 16 -c 1000000 -p sampleA sA_R1.fq sA_R2.fq\n",
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
        DEBUG("Slot available (%d/%d files). Resuming.\n", current_files, max_files);
      }
      return;
    }

    if (!reported_waiting) {
      DEBUG("Max files limit (%d/%d) reached. Waiting...\n", current_files, max_files);
      reported_waiting = true;
    }

    usleep(50000);  // Wait 0.05 seconds
  }
}

// Copy 4 lines from input to output, flushing and refilling buffers as necessary.
// Returns 1 on success, 0 for clean EOF (at the start of a 4-line record), -1 for
// unexpected EOF.
int copy_record(InputState *input, OutputState *output, const char *out_path_for_errors) {
  // input and output must be valid
  if (output->f == NULL) {
    DEBUG("Internal Error: copy_record called with NULL output file pointer.\n");
    exit(1);
  }
  if (input->eof_flag) {
    DEBUG("Called copy_record while at eof.\n");
    exit(1);
  }

  int newline_count = 0;
  while (newline_count < 4) {
    // Refill input buffer if needed
    if (input->pos >= input->valid) {
      // size of 1 since elements are single characters
      input->valid = fread(input->buf, 1, input->buf_size, input->f);
      input->pos = 0;

      if (input->valid == 0) {
        if (feof(input->f)) {
          input->eof_flag = true;
          if (newline_count == 0) return 0;  // clean EOF right at start of record
          fprintf(stderr, "Error: Unexpected EOF in input '%s' after %d lines of record.\n",
                  input->path, newline_count);
          return -1;  // Unexpected EOF mid-record
        } else {      // ferror occurred
          perror_and_exit("fread error in copy_record", input->path);
        }
      }
    }

    // Process characters currently in the input buffer
    while (input->pos < input->valid && newline_count < 4) {
      char c = input->buf[input->pos++];
      output->buf[output->pos++] = c;

      // Flush output buffer if full
      if (output->pos == output->buf_size) {
        // size of 1 for buffer of characters
        if (fwrite(output->buf, 1, output->buf_size, output->f) != output->buf_size) {
          perror_and_exit("fwrite error for full buffer in copy_record",
                          out_path_for_errors);
        }
        output->pos = 0;  // Reset output buffer position
      }

      // Count newline
      if (c == '\n') {
        newline_count++;
      }
    }  // End character processing while
  }  // End while (newline_count < 4)

  return 1;  // Success, 4 lines processed
}
