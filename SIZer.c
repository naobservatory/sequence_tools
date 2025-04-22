#include <errno.h>
#include <getopt.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#define MAX_READ_LEN 1024
#define DEFAULT_TEMP_CHUNK_PREFIX "work_chunk"  // for local uncompressed files

// Struct to hold buffers for reading a FASTQ record
typedef struct {
  char *title;
  size_t title_size;
  char *seq;
  size_t seq_size;
  char *plus;
  size_t plus_size;
  char *quality;
  size_t quality_size;
} FastqRecordBuffers;

// --- Function Prototypes ---

void print_usage_and_exit(const char *program_name, const char *default_prefix);
int reap_children(int *running_children_ptr, int *failed_children_ptr);
void wait_for_slot(int *running_children_ptr, int max_children, int *failed_children_ptr);
void wait_for_all_children(int *running_children_ptr, int *failed_children_ptr);
void check_child_status(pid_t pid, int status, int *failed_children_ptr, const char *context);
pid_t launch_compression_job(const char *local_fname, int chunk_index,
                             const char *s3_prefix, const char *script_path,
                             const char *zstd_threads_str);
FILE *open_fastq_file(const char *path, int *is_pipe);
void close_fastq_file(FILE *f, int is_pipe);
int read_fastq_record(FILE *f, FastqRecordBuffers *buffers, long read_num);
void cleanup_resources(FILE *f1, int p1, FILE *f2, int p2, FastqRecordBuffers *b1,
                       FastqRecordBuffers *b2, FILE *out_chunk_fp);

// --- Main Function ---

/**
 * @brief Main entry point for paired FASTQ -> SIZ.
 *
 * Reads paired-end FASTQ files (R1 and R2), interleaves them, splits them into chunks of a
 * specified number of read pairs, and launches background jobs via an external script to
 * compress each chunk with zstd and upload it to S3. Manages the number of concurrent
 * background jobs.
 *
 * @param argc Argument count.
 * @param argv Argument vector. Expects optional flag -p <prefix> and 7 positional args.
 *  Positional Args (after options):
 *  1: reads_per_chunk (long) - Max read pairs per output chunk.
 *  2: r1_fastq (char*) - Path to R1 FASTQ file (can be .gz).
 *  3: r2_fastq (char*) - Path to R2 FASTQ file (can be .gz).
 *  4: s3_output_prefix (char*) - S3 prefix for output files (e.g., s3://bucket/sample).
 *     Chunks will be named <s3_output_prefix>_chunkNNNNNN.fastq.zst.
 *  5: compress_script_path (char*) - Path to the executable script that handles
 *     compression and upload (e.g., compress_upload.sh).
 *  6: max_concurrent_jobs (int) - Max number of background compression/upload jobs to run
 *     simultaneously.
 *  7: zstd_threads (int) - Number of zstd threads (passed to script).
 *
 * @return int 0 on success (all jobs launched and finished successfully).
 *             1 on error (argument error, file error, fork/exec error, or if any background
 *             jobs failed).
 */
int main(int argc, char **argv) {
  const char *local_prefix = DEFAULT_TEMP_CHUNK_PREFIX;
  int opt;

  // --- Argument Parsing ---
  // Parse command line options using getopt.
  while ((opt = getopt(argc, argv, "p:")) != -1) {
    switch (opt) {
      case 'p':
        local_prefix = optarg;  // Use provided prefix for local temp files
        break;
      default:
        print_usage_and_exit(argv[0], DEFAULT_TEMP_CHUNK_PREFIX);
    }
  }

  // Parse positional arguments.
  // After getopt call above, optind is the index of the first positional arg.
  if (argc - optind != 7) {
    print_usage_and_exit(argv[0], DEFAULT_TEMP_CHUNK_PREFIX);
  }
  long max_reads_per_chunk = strtol(argv[optind], NULL, 10);
  char *r1_path = argv[optind + 1];
  char *r2_path = argv[optind + 2];
  char *s3_prefix = argv[optind + 3];
  char *script_path = argv[optind + 4];
  int max_children = atoi(argv[optind + 5]);
  char *zstd_threads_str = argv[optind + 6];

  // Validate arguments
  if (max_reads_per_chunk <= 0) {
    fprintf(stderr, "Error: reads_per_chunk must be positive\n");
    print_usage_and_exit(argv[0], DEFAULT_TEMP_CHUNK_PREFIX);
  }
  if (max_children <= 0) {
    fprintf(stderr, "Error: max_concurrent_jobs must be positive\n");
    print_usage_and_exit(argv[0], DEFAULT_TEMP_CHUNK_PREFIX);
  }
  if (atoi(zstd_threads_str) <= 0) {
    // We validate here, but pass the string to the script.
    fprintf(stderr, "Error: zstd_threads must be positive\n");
    print_usage_and_exit(argv[0], DEFAULT_TEMP_CHUNK_PREFIX);
  }

  // --- Resource Initialization ---
  FILE *f1 = NULL;               // Input R1 file pointer
  FILE *f2 = NULL;               // Input R2 file pointer
  FILE *out_chunk_fp = NULL;     // Output chunk file pointer
  int f1_pipe = 0, f2_pipe = 0;  // Flags for whether inputs are pipes (popen)

  // Buffers for reading FASTQ records using getline (dynamically allocated)
  FastqRecordBuffers r1_buffers = {NULL, 0, NULL, 0, NULL, 0, NULL, 0};
  FastqRecordBuffers r2_buffers = {NULL, 0, NULL, 0, NULL, 0, NULL, 0};

  // Buffer for temporary local filenames (stack allocated)
  // Size fits prefix length, "_chunk", 6 digits, ".fastq", null terminator.
  int max_local_fname_len = strlen(local_prefix) + 36;
  char local_chunk_fname[max_local_fname_len];

  // Process management variables
  int running_children = 0;
  int failed_children = 0;       // Count of background jobs that failed
  int launched_jobs = 0;         // Total jobs launched
  long current_chunk_reads = 0;  // Reads written to current chunk
  long total_read_pairs = 0;     // Total read pairs processed
  int chunk_index = 0;           // Index for naming chunks

  int final_exit_code = 0;  // Stores the final exit code

  // --- Open Input Files ---
  // open_fastq_file handles .gz via popen and exits fatally on error.
  f1 = open_fastq_file(r1_path, &f1_pipe);
  f2 = open_fastq_file(r2_path, &f2_pipe);

  fprintf(stderr,
          "Starting processing. Max concurrent compression jobs: %d. Temp "
          "prefix: %s\n",
          max_children, local_prefix);

  // --- Main Processing Loop ---
  // Reads pairs, writes to temporary local chunk file, launches background job
  // when chunk is full or EOF is reached.
  while (true) {
    // Read one record pair. read_fastq_record exits fatally on errors.
    int r1_success = read_fastq_record(f1, &r1_buffers, total_read_pairs + 1);
    int r2_success = read_fastq_record(f2, &r2_buffers, total_read_pairs + 1);

    // Determine if we've reached the end of input.
    bool eof_reached = (!r1_success && !r2_success);
    // Check for mismatched reads (one file ends before the other).
    if (!eof_reached && (r1_success != r2_success)) {
      fprintf(stderr,
              "Error: Unequal number of reads in R1 and R2 files (around read "
              "pair #%ld).\n",
              total_read_pairs + 1);
      final_exit_code = 1;
      break;  // Exit loop to proceed to cleanup and final wait
    }

    // If not EOF, process the read pair.
    if (!eof_reached) {
      // Open new output chunk file if it's not already open.
      if (!out_chunk_fp) {
        snprintf(local_chunk_fname, max_local_fname_len, "%s_chunk%06d.fastq", local_prefix,
                 chunk_index);
        out_chunk_fp = fopen(local_chunk_fname, "w");
        if (!out_chunk_fp) {
          fprintf(stderr, "Error creating output chunk file %s: %s\n", local_chunk_fname,
                  strerror(errno));
          final_exit_code = 1;
          break;  // Exit loop for cleanup
        }
        fprintf(stderr, "Opened chunk file %s\n", local_chunk_fname);
      }

      // Write interleaved records to current chunk file.
      fprintf(out_chunk_fp, "%s%s%s%s", r1_buffers.title, r1_buffers.seq, r1_buffers.plus,
              r1_buffers.quality);
      fprintf(out_chunk_fp, "%s%s%s%s", r2_buffers.title, r2_buffers.seq, r2_buffers.plus,
              r2_buffers.quality);

      total_read_pairs++;
      current_chunk_reads++;
    }

    // Check if current chunk is full OR if we reached EOF and have data in the
    // current chunk.
    bool chunk_ready_to_launch =
        (current_chunk_reads > 0 &&
         (current_chunk_reads >= max_reads_per_chunk || eof_reached));

    if (chunk_ready_to_launch) {
      if (fclose(out_chunk_fp) != 0) {
        perror("Warning: fclose failed for chunk file");
        // Continue to launch job, but data might be incomplete.
        // Consider marking as failed? For now, just warn.
      }
      out_chunk_fp = NULL;  // Mark as closed
      fprintf(stderr, "Closed chunk file %s (%ld reads)\n", local_chunk_fname,
              current_chunk_reads);

      // Wait for an available slot if the max number of background jobs are running.
      wait_for_slot(&running_children, max_children, &failed_children);

      // Launch the background compression/upload job.
      pid_t pid = launch_compression_job(local_chunk_fname, chunk_index, s3_prefix,
                                         script_path, zstd_threads_str);
      if (pid < 0) {  // Fork or other error in launch_compression_job
        fprintf(stderr, "Error: Failed to launch job for chunk %d.\n", chunk_index);
        // Don't try to launch more jobs. Proceed to cleanup.
        final_exit_code = 1;
        break;  // Exit loop
      } else {
        running_children++;
        launched_jobs++;
        fprintf(stderr, "Parent: Launched job for chunk %d (PID %d) (%d running)\n",
                chunk_index, pid, running_children);
      }

      // Reset for next chunk
      current_chunk_reads = 0;
      chunk_index++;
    }

    // If we reached EOF, break the main loop after potentially launching the last job.
    if (eof_reached) {
      break;
    }
  }  // End while(true) main processing loop

  // --- Cleanup and Final Wait ---
  // Free getline buffers, close input files, handle any lingering open chunk file.
  cleanup_resources(f1, f1_pipe, f2, f2_pipe, &r1_buffers, &r2_buffers, out_chunk_fp);
  f1 = f2 = out_chunk_fp = NULL;  // Avoid accidental reuse
  wait_for_all_children(&running_children, &failed_children);
  fprintf(stderr, "Processed %ld read pairs. Launched %d compression/upload jobs.\n",
          total_read_pairs, launched_jobs);

  // Set final exit code based on earlier errors OR failed background jobs.
  if (failed_children > 0) {
    fprintf(stderr, "Error: %d background job(s) failed.\n", failed_children);
    final_exit_code = 1;  // Ensure failure is indicated
  }
  fprintf(stderr, "All background jobs accounted for. Exiting with code %d.\n",
          final_exit_code);
  return final_exit_code;
}

// --- Helper Function Definitions ---

/**
 * @brief Prints usage instructions and exits.
 *
 * @param program_name The name of the program (argv[0]).
 * @param default_prefix The default temporary file prefix.
 */
void print_usage_and_exit(const char *program_name, const char *default_prefix) {
  fprintf(stderr,
          "Usage: %s [-p <local_prefix>] <reads_per_chunk> <r1_fastq> "
          "<r2_fastq> <s3_output_prefix> <compress_script_path> "
          "<max_concurrent_jobs> <zstd_threads>\n",
          program_name);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -p <local_prefix>  Prefix for temporary chunk files (default: %s)\n",
          default_prefix);
  fprintf(stderr, "Arguments:\n");
  fprintf(stderr, "  reads_per_chunk      Max read pairs per output chunk\n");
  fprintf(stderr, "  r1_fastq             Path to R1 FASTQ file (can be .gz)\n");
  fprintf(stderr, "  r2_fastq             Path to R2 FASTQ file (can be .gz)\n");
  fprintf(stderr, "  s3_output_prefix     S3 prefix for output files\n");
  fprintf(stderr, "  compress_script_path Path to compression/upload script\n");
  fprintf(stderr, "  max_concurrent_jobs  Max simultaneous background jobs\n");
  fprintf(stderr, "  zstd_threads         Number of zstd threads for script\n");
  fprintf(stderr, "\nExample:\n");
  fprintf(stderr,
          "  %s -p work/sampleA 1000000 r1.fq.gz r2.fq.gz "
          "s3://mybucket/siz/sampleA scripts/compress_upload.sh 16 3\n",
          program_name);
  exit(1);
}

/**
 * @brief Cleans up allocated resources like file pointers and buffers.
 * Intended to be called before exiting main.
 *
 * @param f1 R1 input file pointer.
 * @param p1 R1 pipe flag.
 * @param f2 R2 input file pointer.
 * @param p2 R2 pipe flag.
 * @param b1 R1 buffer struct pointer.
 * @param b2 R2 buffer struct pointer.
 * @param out_chunk_fp Output chunk file pointer (may be NULL or open).
 */
void cleanup_resources(FILE *f1, int p1, FILE *f2, int p2, FastqRecordBuffers *b1,
                       FastqRecordBuffers *b2, FILE *out_chunk_fp) {
  // Close input files (safe to call even if NULL)
  close_fastq_file(f1, p1);
  close_fastq_file(f2, p2);

  // Free all line buffers allocated by getline (safe to call on NULL)
  free(b1->title);
  free(b1->seq);
  free(b1->plus);
  free(b1->quality);
  free(b2->title);
  free(b2->seq);
  free(b2->plus);
  free(b2->quality);

  // Ensure any unclosed output chunk file is handled (e.g., if loop exited abnormally)
  // We don't launch a compress-upload job for this potentially-corrupt file.
  if (out_chunk_fp != NULL) {
    fprintf(stderr,
            "Warning: Output chunk file was open during cleanup. Closing "
            "now (data may be incomplete).\n");
    if (fclose(out_chunk_fp) != 0) {
      perror("fclose failed during cleanup");
    }
  }
}

/**
 * @brief Helper function to check and report child process status.
 * Updates the failed children count if the process didn't exit normally with status 0.
 *
 * @param pid The process ID of the child that was reaped.
 * @param status The status value returned by wait/waitpid.
 * @param failed_children_ptr Pointer to the failed children count to update.
 * @param context A string describing the context for logging purposes.
 */
void check_child_status(pid_t pid, int status, int *failed_children_ptr,
                               const char *context) {
  if (WIFEXITED(status) && WEXITSTATUS(status) == 0) {
    // Child exited normally with status 0 (success)
    fprintf(stderr, "%s: Child %d finished successfully.\n", context, pid);
  } else {
    // Child either didn't exit normally or exited with non-zero status
    if (WIFEXITED(status)) {
      fprintf(stderr, "Warning: %s: Child %d exited with status %d.\n", context, pid,
              WEXITSTATUS(status));
    } else if (WIFSIGNALED(status)) {
      fprintf(stderr, "Warning: %s: Child %d killed by signal %d.\n", context, pid,
              WTERMSIG(status));
    } else {
      fprintf(stderr, "Warning: %s: Child %d finished abnormally.\n", context, pid);
    }
    (*failed_children_ptr)++;
  }
}

/**
 * @brief Non-blockingly checks for and reaps any terminated child processes.
 * Updates the count of running children and failed children. This is used
 * periodically to check for finished jobs without blocking.
 *
 * @param running_children_ptr Pointer to the running children count. Decremented
 * on reap.
 * @param failed_children_ptr Pointer to the failed children count. Incremented
 * on failure.
 * @return int The number of children successfully reaped in this call.
 */
int reap_children(int *running_children_ptr, int *failed_children_ptr) {
  int reaped_count = 0;
  pid_t finished_pid;
  int status;

  // Loop as long as waitpid finds finished children without blocking (WNOHANG)
  while ((finished_pid = waitpid(-1, &status, WNOHANG)) > 0) {
    (*running_children_ptr)--;
    reaped_count++;
    check_child_status(finished_pid, status, failed_children_ptr, "Non-blocking reap");
  }

  // Handle errors from waitpid itself, ignoring ECHILD (no children left)
  if (finished_pid == -1 && errno != ECHILD) {
    perror("Warning: waitpid error in reap_children");
  }
  return reaped_count;
}

/**
 * @brief Waits until a slot is available for a new background job. It first tries
 * non-blocking reaping, then blocks if the job limit is still reached.
 *
 * @param running_children_ptr Pointer to the running children count.
 * @param max_children The maximum allowed running children.
 * @param failed_children_ptr Pointer to the failed children count.
 */
void wait_for_slot(int *running_children_ptr, int max_children, int *failed_children_ptr) {
  // First, try to reap any finished children without blocking.
  reap_children(running_children_ptr, failed_children_ptr);

  // If we're still at or above the limit, we need to wait blockingly.
  while (*running_children_ptr >= max_children) {
    fprintf(stderr, "Reached max concurrent jobs (%d). Waiting for a slot...\n",
            max_children);
    pid_t finished_pid;
    int status;

    // Wait blockingly for ANY child process to change state.
    // Retry if interrupted by a signal (EINTR).
    do {
      finished_pid = wait(&status);
    } while (finished_pid == -1 && errno == EINTR);

    if (finished_pid < 0) {
      // Handle unexpected errors from wait().
      if (errno == ECHILD) {
        // This shouldn't happen if running_children_ptr >= max_children > 0,
        // but handle defensively.
        fprintf(stderr,
                "Warning: wait() returned ECHILD unexpectedly in "
                "wait_for_slot. Resetting count.\n");
        *running_children_ptr = 0;  // Reset count to avoid infinite loop
        break;                      // Exit the waiting loop
      } else {
        // For other errors (EPERM, etc.), it's likely fatal.
        perror("Fatal: wait error while waiting for slot");
        exit(1);
      }
    } else {
      // Successfully waited for a child.
      (*running_children_ptr)--;

      // Use the common helper to check and report status
      check_child_status(finished_pid, status, failed_children_ptr, "Slot wait");
    }
  }
}

/**
 * @brief Waits blockingly for all remaining child processes to complete.
 * Called at the end of the program before exiting.
 *
 * @param running_children_ptr Pointer to the count of running children. Will be set to 0
 * when complete.
 * @param failed_children_ptr Pointer to the count of failed children. Will be incremented
 * for each child that fails.
 */
void wait_for_all_children(int *running_children_ptr, int *failed_children_ptr) {
  fprintf(stderr,
          "Cleanup phase: Waiting for %d remaining background jobs to "
          "finish...\n",
          *running_children_ptr);

  // Loop while we know there are children running.
  while (*running_children_ptr > 0) {
    pid_t finished_pid;
    int status;

    // Wait blockingly for ANY child process to change state.
    // Retry if interrupted by a signal (EINTR).
    do {
      finished_pid = wait(&status);
    } while (finished_pid == -1 && errno == EINTR);

    if (finished_pid < 0) {
      // Handle unexpected errors from wait().
      if (errno == ECHILD) {
        // We thought children were running, but wait() says no. Log it.
        fprintf(stderr,
                "Warning: wait() returned ECHILD during final cleanup, but "
                "expected %d children. Resetting count.\n",
                *running_children_ptr);
      } else {
        perror("Warning: wait error during final cleanup");
        // Consider this scenario as indicating a failure.
        (*failed_children_ptr)++;
      }
      *running_children_ptr = 0;  // Reset the count to avoid infinite loop
    } else {
      // Successfully waited for a child.
      (*running_children_ptr)--;

      // Use the common helper to check and report status
      check_child_status(finished_pid, status, failed_children_ptr, "Final wait");
    }
  }
}

/**
 * @brief Creates a child process to execute the compression/upload script.
 * Handles formatting the S3 path argument for the script.
 * Child process exits(1) if exec fails or path formatting fails.
 * Parent process returns PID of child or -1 on fork error.
 *
 * @param local_fname Path to the local uncompressed chunk file.
 * @param chunk_index The numerical index of this chunk.
 * @param s3_prefix S3 prefix for the final output object.
 * @param script_path Path to the executable compression/upload script.
 * @param zstd_threads_str Number of zstd threads (as string) for the script.
 * @return pid_t Child PID on success, or -1 if fork() failed in parent.
 */
pid_t launch_compression_job(const char *local_fname, int chunk_index,
                             const char *s3_prefix, const char *script_path,
                             const char *zstd_threads_str) {
  pid_t pid = fork();
  if (pid < 0) {
    // Error occurred in the parent process during fork.
    perror("fork failed in launch_compression_job");
    return -1;
  }

  if (pid == 0) {  // Child process
    // Construct the target S3 path string.
    // Buffer must fit: prefix + "_chunk" + 6 digits + ".fastq.zst" + null
    int s3_fname_buf_size = strlen(s3_prefix) + 36;
    char s3_fname[s3_fname_buf_size];
    int written = snprintf(s3_fname, s3_fname_buf_size, "%s_chunk%06d.fastq.zst", s3_prefix,
                           chunk_index);

    // Check for snprintf errors or truncation
    if (written < 0 || written >= s3_fname_buf_size) {
      fprintf(stderr,
              "Child %d: Error formatting S3 path (written=%d, size=%d). "
              "Prefix or index too long?\n",
              getpid(), written, s3_fname_buf_size);
      exit(1);  // Exit child cleanly on formatting error
    }

    // Prepare arguments for execlp. The first argument is the script path,
    // the following are argv[0], argv[1], ... for the script, terminated by NULL.
    fprintf(stderr, "Child %d: Launching script for %s -> %s\n", getpid(), local_fname,
            s3_fname);
    execlp(script_path, script_path, local_fname, s3_fname, zstd_threads_str, (char *)NULL);

    // If execlp returns, it means execution failed.
    fprintf(stderr, "Child %d: Failed to exec script '%s': %s\n", getpid(), script_path,
            strerror(errno));
    exit(1);  // Exit child indicating failure
  } else {    // Parent process
    // Successfully launched child.
    return pid;  // Return child PID to the main process
  }
}

/**
 * @brief Opens a FASTQ file, handling .gz compression via popen/zcat if needed.
 * Exits fatally if the file cannot be opened or the zcat command fails.
 *
 * @param path Path to the FASTQ file (.gz or plain).
 * @param is_pipe Output pointer; set to true if popen was used, false otherwise.
 * @return FILE* A file pointer (either from fopen or popen) ready for reading.
 */
FILE *open_fastq_file(const char *path, int *is_pipe) {
  FILE *f;
  // Check if the path ends with ".gz"
  const char *extension = strrchr(path, '.');
  if (extension && extension != path && strcmp(extension, ".gz") == 0) {
    // Construct the "zcat <path>" command.
    // Allocate a buffer large enough for "zcat " + path + null terminator
    char cmd[2048];  // should suffice for all path lengths encountered in practice
    int written = snprintf(cmd, sizeof(cmd), "zcat %s", path);
    // Check for formatting errors or truncation.
    if (written < 0 || (size_t)written >= sizeof(cmd)) {
      fprintf(stderr, "Error: Failed to format zcat command (path likely too long): %s\n",
              path);
      exit(1);
    }
    // Open a pipe to the zcat command for reading.
    f = popen(cmd, "r");
    if (!f) {
      // popen can fail if fork/pipe fails or if the command doesn't exist.
      fprintf(stderr, "Error: popen failed for command '%s': %s\n", cmd, strerror(errno));
      exit(1);
    }
    *is_pipe = true;
  } else {
    // Open as a regular file.
    f = fopen(path, "r");
    if (!f) {
      fprintf(stderr, "Error opening input file %s: %s\n", path, strerror(errno));
      exit(1);
    }
    *is_pipe = false;
  }
  return f;
}

/**
 * @brief Closes a file pointer previously opened by open_fastq_file.
 *
 * @param f The file pointer to close. Can be NULL (no-op).
 * @param is_pipe Flag indicating if the file pointer originated from popen.
 */
void close_fastq_file(FILE *f, int is_pipe) {
  if (!f) return;  // Do nothing if the pointer is NULL

  if (is_pipe) {
    int status = pclose(f);  // Close pipe and get exit status of command
    if (status == -1) {
      // pclose() failed itself.
      perror("Warning: pclose failed on input stream");
    } else if (WIFEXITED(status) && WEXITSTATUS(status) != 0) {
      // Command (e.g., zcat) exited with non-zero status.
      fprintf(stderr, "Warning: Input pipe command finished with exit status %d\n",
              WEXITSTATUS(status));
    } else if (WIFSIGNALED(status)) {
      // Command was terminated by a signal.
      fprintf(stderr, "Warning: Input pipe command killed by signal %d\n",
              WTERMSIG(status));
    }
    // We generally don't treat input stream close errors as fatal for the
    // overall process.
  } else {
    if (fclose(f) != 0) {
      // fclose failed for a regular file.
      perror("Warning: fclose failed on input stream");
    }
  }
}

/**
 * @brief Reads one complete FASTQ record (4 lines) from the file stream.
 * Uses getline for dynamic buffer allocation. Performs basic validation.
 * Exits fatally on read errors or format corruption.
 * Does not support hard wrapping.
 *
 * @param f Input file stream (from fopen or popen).
 * @param buffers Pointer to FastqRecordBuffers struct for storing lines.
 * @param read_num The 1-based index of the read pair being processed (for error
 * messages).
 * @return int 1 on successful read, 0 on clean EOF before starting the record.
 */
int read_fastq_record(FILE *f, FastqRecordBuffers *buffers, long read_num) {
  // Use getline to read each of the four lines. getline handles buffer resizing.
  ssize_t title_len = getline(&(buffers->title), &(buffers->title_size), f);
  // Check for EOF or error immediately after the first line read.
  if (title_len == -1) {
    if (feof(f)) {
      return 0;
    }  // Clean EOF before record started
    else {  // Read error
      perror("getline error reading FASTQ title");
      fprintf(stderr, "Error occurred reading title for read #%ld.\n", read_num);
      exit(1);
    }
  }

  ssize_t seq_len = getline(&(buffers->seq), &(buffers->seq_size), f);
  ssize_t plus_len = getline(&(buffers->plus), &(buffers->plus_size), f);
  ssize_t quality_len = getline(&(buffers->quality), &(buffers->quality_size), f);

  // Check if any subsequent getline failed.
  if (seq_len == -1 || plus_len == -1 || quality_len == -1) {
    if (feof(f)) {
      // Reached EOF unexpectedly after already reading part of a record.
      fprintf(stderr, "Error: Incomplete FASTQ record at end of file (Read #%ld).\n",
              read_num);
    } else {
      // Read error occurred mid-record.
      perror("getline error reading FASTQ record");
      fprintf(stderr, "Error occurred reading data for read #%ld.\n", read_num);
    }
    exit(1);  // Treat incomplete records or read errors as fatal
  }

  // Basic FASTQ format validation
  // Check start characters, sequence/quality length match, and reasonable line length.
  // Note: seq_len includes the newline, so it should equal quality_len.
  if ((buffers->title)[0] != '@' || (buffers->plus)[0] != '+' || seq_len != quality_len ||
      title_len > MAX_READ_LEN || seq_len > MAX_READ_LEN || plus_len > MAX_READ_LEN ||
      quality_len > MAX_READ_LEN) {
    fprintf(stderr, "Error: Corrupt FASTQ record detected processing read #%ld\n",
            read_num);
    fprintf(stderr,
            "       Record starts: Title='%c' Plus='%c'. Lengths: Seq=%zd "
            "Qual=%zd\n",
            (buffers->title)[0], (buffers->plus)[0], seq_len, quality_len);
    exit(1);  // Treat corruption as fatal
  }

  return 1;  // Successfully read a full, valid record
}
