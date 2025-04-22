#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>  // For pid_t
#include <sys/wait.h>
#include <unistd.h>

#define MAX_READ_LEN 1024
#define TEMP_CHUNK_PREFIX \
  "work_chunk"  // Prefix for temporary local uncompressed files

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

int reap_children(int *running_children_ptr, int *failed_children_ptr);
void wait_for_slot(int *running_children_ptr, int max_children,
                   int *failed_children_ptr);
pid_t launch_compression_job(const char *local_fname, int chunk_index,
                             const char *s3_prefix, const char *script_path,
                             const char *zstd_threads_str);
FILE *open_fastq_file(const char *path, int *is_pipe);
void close_fastq_file(FILE *f, int is_pipe);
int read_fastq_record(FILE *f, FastqRecordBuffers *buffers, long read_num);

/**
 * @brief Main entry point for the FASTQ chunking and compression program.
 *
 * Reads paired-end FASTQ files (R1 and R2), interleaves them, splits them
 * into chunks of a specified number of read pairs, and launches background
 * jobs via an external script to compress each chunk with zstd and upload it to
 * S3. Manages the number of concurrent background jobs.
 *
 * @param argc Argument count. Expected 7 arguments + program name.
 * @param argv Argument vector:
 * argv[1]: reads_per_chunk (long): Max read pairs per output chunk.
 * argv[2]: r1_fastq (char*): Path to R1 FASTQ file (can be .gz).
 * argv[3]: r2_fastq (char*): Path to R2 FASTQ file (can be .gz).
 * argv[4]: s3_output_prefix (char*): S3 prefix for output files
 * (e.g., s3://bucket/path/sample). Chunks will be named
 * <s3_output_prefix>_chunkNNNNNN.fastq.zst.
 * argv[5]: compress_script_path (char*): Path to the executable script
 * that handles compression and upload (e.g., compress_upload.sh).
 * argv[6]: max_concurrent_jobs (int): Max number of background
 * compression/upload jobs to run simultaneously. argv[7]: zstd_threads (char*):
 * Number of zstd threads (passed to script).
 *
 * @return int 0 on success (all jobs launched and finished successfully).
 * 1 on error (argument error, file error, fork/exec error, or if any
 * background jobs failed).
 */
int main(int argc, char **argv) {
  // Use exit(1) directly for argument parsing errors, as no resources allocated
  // yet.
  if (argc != 8) {
    fprintf(
        stderr,
        "Usage: %s <reads_per_chunk> <r1_fastq> <r2_fastq> <s3_output_prefix> "
        "<compress_script_path> <max_concurrent_jobs> <zstd_threads>\n",
        argv[0]);
    fprintf(stderr,
            "  Example: %s 1000000 r1.fq.gz r2.fq.gz s3://mybucket/siz/sample1 "
            "scripts/compress_upload.sh 16 3\n",
            argv[0]);
    fprintf(stderr, "  Generates local temp files named %s_chunkNNNNNN.fastq\n",
            TEMP_CHUNK_PREFIX);
    fprintf(stderr,
            "  Launches <compress_script_path> for each chunk, targeting "
            "<s3_output_prefix>_chunkNNNNNN.fastq.zst\n");
    exit(1);
  }

  // Parse arguments
  long max_reads_per_chunk = strtol(argv[1], NULL, 10);
  char *r1_path = argv[2];
  char *r2_path = argv[3];
  char *s3_prefix = argv[4];
  char *script_path = argv[5];
  int max_children = atoi(argv[6]);
  char *zstd_threads_str = argv[7];

  // Validate arguments
  if (max_reads_per_chunk <= 0) {
    fprintf(stderr, "Error: reads_per_chunk must be positive\n");
    exit(1);
  }
  if (max_children <= 0) {
    fprintf(stderr, "Error: max_concurrent_jobs must be positive\n");
    exit(1);
  }
  if (atoi(zstd_threads_str) <= 0) {
    fprintf(stderr, "Error: zstd_threads must be positive\n");
    exit(1);
  }

  // --- Resource Allocation Starts Here ---
  int exit_code = 0;  // Final exit code for main, assume success initially
  FILE *f1 = NULL;    // Initialize file pointers
  FILE *f2 = NULL;
  FILE *out_chunk_fp = NULL;
  int f1_pipe = 0, f2_pipe = 0;

  // Buffers for reading FASTQ records - initialized safely by struct definition
  FastqRecordBuffers r1_buffers = {NULL, 0, NULL, 0, NULL, 0, NULL, 0};
  FastqRecordBuffers r2_buffers = {NULL, 0, NULL, 0, NULL, 0, NULL, 0};

  // Buffer for temporary local filenames (stack allocated, contents undefined
  // until snprintf)
  int max_local_fname_len = strlen(TEMP_CHUNK_PREFIX) + 36;
  char local_chunk_fname[max_local_fname_len];

  // Process management variables
  int running_children = 0;
  int failed_children = 0;
  int launched_jobs = 0;
  long current_chunk_reads = 0;
  long total_read_pairs = 0;
  int chunk_index = 0;

  // Open input files (exit(1) here is acceptable as child processes/buffers not
  // active yet)
  f1 = open_fastq_file(r1_path, &f1_pipe);
  f2 = open_fastq_file(
      r2_path,
      &f2_pipe);  // If this fails after f1 opens, cleanup_exit handles f1

  fprintf(stderr, "Starting processing. Max concurrent compression jobs: %d\n",
          max_children);

  // Main processing loop: read pairs, write to chunk, launch job when chunk is
  // full
  while (1) {
    // Read one record pair
    // read_fastq_record exits fatally on read error or corrupt data
    int r1_success = read_fastq_record(f1, &r1_buffers, total_read_pairs);
    int r2_success = read_fastq_record(f2, &r2_buffers, total_read_pairs);

    // Check for EOF or mismatched reads
    if (!r1_success && !r2_success) {
      break;
    }  // Normal EOF reached on both files
    if (r1_success != r2_success) {
      fprintf(stderr,
              "Error: Unequal number of reads in R1 and R2 files (around read "
              "pair #%ld).\n",
              total_read_pairs + 1);
      exit_code = 1;      // Mark failure
      goto cleanup_exit;  // Use goto for central cleanup on error
    }

    // Open new output chunk file if necessary
    if (!out_chunk_fp) {
      snprintf(local_chunk_fname, max_local_fname_len, "%s_chunk%06d.fastq",
               TEMP_CHUNK_PREFIX, chunk_index);
      out_chunk_fp = fopen(local_chunk_fname, "w");
      if (!out_chunk_fp) {
        fprintf(stderr, "Error creating output chunk file %s: %s\n",
                local_chunk_fname, strerror(errno));
        exit_code = 1;      // Mark failure
        goto cleanup_exit;  // Use goto as file pointers might be open
      }
      fprintf(stderr, "Opened chunk file %s\n", local_chunk_fname);
    }

    // Write interleaved records to current chunk file
    fprintf(out_chunk_fp, "%s%s%s%s", r1_buffers.title, r1_buffers.seq,
            r1_buffers.plus, r1_buffers.quality);
    fprintf(out_chunk_fp, "%s%s%s%s", r2_buffers.title, r2_buffers.seq,
            r2_buffers.plus, r2_buffers.quality);

    total_read_pairs++;
    current_chunk_reads++;

    // Check if current chunk is full
    if (current_chunk_reads >= max_reads_per_chunk) {
      if (fclose(out_chunk_fp) != 0) {
        perror("fclose failed for chunk file");
        // Decide if this is fatal? For now, log and continue to launch job.
      }
      out_chunk_fp = NULL;  // Mark as closed
      fprintf(stderr, "Closed chunk file %s (%ld reads)\n", local_chunk_fname,
              current_chunk_reads);

      // Wait for an available slot if max concurrent jobs are running
      wait_for_slot(&running_children, max_children, &failed_children);

      // Launch the background compression/upload job
      pid_t pid =
          launch_compression_job(local_chunk_fname, chunk_index, s3_prefix,
                                 script_path, zstd_threads_str);
      if (pid < 0) {  // Fork failed in helper
        fprintf(stderr, "Error: Failed to launch job for chunk %d.\n",
                chunk_index);
        exit_code = 1;      // Mark failure
        goto cleanup_exit;  // Use goto as children might be running
      } else {
        running_children++;
        launched_jobs++;
        fprintf(stderr,
                "Parent: Launched job for chunk %d (PID %d) (%d running)\n",
                chunk_index, pid, running_children);
      }

      // Reset for next chunk
      current_chunk_reads = 0;
      chunk_index++;
    }  // End if chunk is full
  }  // End while(1) main processing loop

  // Handle the final partial chunk if it contains data
  if (out_chunk_fp) {
    if (fclose(out_chunk_fp) != 0) {
      perror("fclose failed for final chunk file");
    }
    out_chunk_fp = NULL;  // Ensure it's marked closed
    fprintf(stderr, "Closed final chunk file %s (%ld reads)\n",
            local_chunk_fname, current_chunk_reads);

    // Wait for slot and launch job for the final chunk
    wait_for_slot(&running_children, max_children, &failed_children);
    pid_t pid =
        launch_compression_job(local_chunk_fname, chunk_index, s3_prefix,
                               script_path, zstd_threads_str);
    if (pid < 0) {
      fprintf(stderr, "Error: Failed to launch job for final chunk %d.\n",
              chunk_index);
      exit_code = 1;  // Mark failure
                      // Proceed to cleanup as file is closed and error logged
    } else {
      running_children++;
      launched_jobs++;
      fprintf(stderr,
              "Parent: Launched job for final chunk %d (PID %d) (%d running)\n",
              chunk_index, pid, running_children);
    }
  }  // End handling final chunk

cleanup_exit:;  // Label for centralized cleanup before exiting (used by goto on
                // error)

  // --- Resource Cleanup Starts Here ---

  // Close input files (safe to call even if NULL or already closed)
  close_fastq_file(f1, f1_pipe);
  close_fastq_file(f2, f2_pipe);
  f1 = f2 = NULL;  // Prevent potential use after close/free

  // Free all line buffers allocated by getline (safe to call on NULL)
  free(r1_buffers.title);
  free(r1_buffers.seq);
  free(r1_buffers.plus);
  free(r1_buffers.quality);
  free(r2_buffers.title);
  free(r2_buffers.seq);
  free(r2_buffers.plus);
  free(r2_buffers.quality);

  // Ensure any unclosed output chunk file is handled (should be NULL if loop
  // exited normally)
  if (out_chunk_fp != NULL) {
    fprintf(stderr,
            "Warning: Output chunk file was not closed before cleanup. Closing "
            "now.\n");
    if (fclose(out_chunk_fp) != 0) {
      perror("fclose failed during cleanup");
      // Don't launch job, data might be incomplete. Mark failure.
      exit_code = 1;
    }
  }

  // Wait for all remaining child processes to complete
  fprintf(
      stderr,
      "Cleanup phase: Waiting for %d remaining background jobs to finish...\n",
      running_children);
  while (running_children > 0) {
    pid_t finished_pid;
    int status;
    // Wait blockingly, retrying if interrupted by a signal
    do {
      finished_pid = wait(&status);
    } while (finished_pid == -1 && errno == EINTR);

    if (finished_pid < 0) {
      // Handle unexpected errors from wait()
      if (errno == ECHILD && running_children > 0) {
        fprintf(stderr,
                "Warning: wait() returned ECHILD during final cleanup, but "
                "expected %d children. Resetting count.\n",
                running_children);
      } else {
        perror("wait error during final cleanup");
      }
      // Mark overall process as failed if unexpected error occurred
      exit_code = 1;
      running_children =
          0;  // Avoid potential infinite loop on persistent error
    } else {
      // Successfully waited for a child
      running_children--;
      fprintf(stderr, "Waited for child %d (%d remaining).\n", finished_pid,
              running_children);
      // Check status and increment failed_children counter if job failed
      if (WIFEXITED(status)) {
        if (WEXITSTATUS(status) != 0) {
          fprintf(stderr, "Warning: Background job %d exited with status %d.\n",
                  finished_pid, WEXITSTATUS(status));
          failed_children++;
        }
      } else if (WIFSIGNALED(status)) {
        fprintf(stderr, "Warning: Background job %d killed by signal %d.\n",
                finished_pid, WTERMSIG(status));
        failed_children++;
      } else {
        fprintf(stderr,
                "Warning: Background job %d finished with unexpected status "
                "code %d.\n",
                finished_pid, status);
        failed_children++;
      }
    }
  }

  fprintf(stderr,
          "Processed %ld read pairs. Launched %d compression/upload jobs.\n",
          total_read_pairs, launched_jobs);

  // Set final exit code based on whether any children failed *or* if an error
  // occurred earlier
  if (failed_children > 0) {
    fprintf(stderr, "Error: %d background job(s) failed. See warnings above.\n",
            failed_children);
    exit_code = 1;  // Ensure failure is indicated
  }

  fprintf(stderr, "All background jobs accounted for. Exiting with code %d.\n",
          exit_code);
  return exit_code;
}

// --- Helper Function Definitions ---

/**
 * @brief Non-blockingly checks for and reaps any terminated child processes.
 * Updates the count of running children and failed children based on exit
 * status.
 *
 * @param running_children_ptr Pointer to the integer tracking the number of
 * running children. This value is decremented when a child is reaped.
 * @param failed_children_ptr Pointer to the integer tracking the number of
 * children that terminated with non-zero status or via signal. This value is
 * incremented on child failure.
 * @return int The number of children successfully reaped in this call.
 */
int reap_children(int *running_children_ptr, int *failed_children_ptr) {
  int reaped_count = 0;
  pid_t finished_pid;
  int status;
  while ((finished_pid = waitpid(-1, &status, WNOHANG)) > 0) {
    (*running_children_ptr)--;
    reaped_count++;
    fprintf(stderr, "Reaped child %d (%d still running).\n", finished_pid,
            *running_children_ptr);
    if (WIFEXITED(status)) {
      if (WEXITSTATUS(status) != 0) {
        fprintf(stderr, "Warning: Reaped child %d exited with status %d.\n",
                finished_pid, WEXITSTATUS(status));
        (*failed_children_ptr)++;
      }
    } else if (WIFSIGNALED(status)) {
      fprintf(stderr, "Warning: Reaped child %d killed by signal %d.\n",
              finished_pid, WTERMSIG(status));
      (*failed_children_ptr)++;
    } else {
      fprintf(stderr,
              "Warning: Reaped child %d finished with unexpected status.\n",
              finished_pid);
      (*failed_children_ptr)++;
    }
  }
  if (finished_pid == -1 && errno != ECHILD) {
    perror("waitpid error in reap_children");
  }
  return reaped_count;
}

/**
 * @brief Waits (blockingly) until the number of running children is less than
 * max_children. Handles EINTR errors by retrying wait(). Updates running/failed
 * counts as children finish. Exits fatally on persistent wait() errors other
 * than EINTR/ECHILD.
 *
 * @param running_children_ptr Pointer to the running children count.
 * @param max_children The maximum allowed running children.
 * @param failed_children_ptr Pointer to the failed children count.
 */
void wait_for_slot(int *running_children_ptr, int max_children,
                   int *failed_children_ptr) {
  reap_children(running_children_ptr,
                failed_children_ptr);  // Try non-blocking first
  while (*running_children_ptr >= max_children) {
    fprintf(stderr, "Reached max concurrent jobs (%d). Waiting for a slot...\n",
            max_children);
    pid_t finished_pid;
    int status;
    do {
      finished_pid = wait(&status);
    } while (finished_pid == -1 && errno == EINTR);  // Retry wait on EINTR
    if (finished_pid < 0) {
      if (errno == ECHILD) {
        fprintf(stderr,
                "Warning: wait() returned ECHILD in wait_for_slot, but "
                "running_children=%d. Resetting count.\n",
                *running_children_ptr);
        *running_children_ptr = 0;
        break;
      } else {
        perror("wait error while waiting for slot");
        exit(1);
      }  // Fail fast on other errors
    } else {
      (*running_children_ptr)--;
      fprintf(stderr, "Slot freed up by child %d (%d still running).\n",
              finished_pid, *running_children_ptr);
      if (WIFEXITED(status)) {
        if (WEXITSTATUS(status) != 0) {
          fprintf(
              stderr,
              "Warning: Slot freed by child %d which exited with status %d.\n",
              finished_pid, WEXITSTATUS(status));
          (*failed_children_ptr)++;
        }
      } else if (WIFSIGNALED(status)) {
        fprintf(
            stderr,
            "Warning: Slot freed by child %d which was killed by signal %d.\n",
            finished_pid, WTERMSIG(status));
        (*failed_children_ptr)++;
      } else {
        fprintf(stderr,
                "Warning: Slot freed by child %d which finished with "
                "unexpected status.\n",
                finished_pid);
        (*failed_children_ptr)++;
      }
    }
  }
}

/**
 * @brief Creates a child process to execute the compression/upload script.
 * Checks snprintf success.
 *
 * @param local_fname Path to the local uncompressed chunk file.
 * @param chunk_index The numerical index of this chunk (used for naming).
 * @param s3_prefix S3 prefix for the final output object.
 * @param script_path Path to the executable compression/upload script.
 * @param zstd_threads_str Number of zstd threads (as a string) to pass to the
 * script.
 * @return pid_t The PID of the launched child process on success, or -1 if
 * fork() failed. Note: Does not wait for the child. The child process exits(1)
 * if exec or snprintf fails.
 */
pid_t launch_compression_job(const char *local_fname, int chunk_index,
                             const char *s3_prefix, const char *script_path,
                             const char *zstd_threads_str) {
  pid_t pid = fork();
  if (pid < 0) {
    perror("fork failed in launch_compression_job");
    return -1;
  }
  if (pid == 0) {                                    // Child process
    int s3_fname_buf_size = strlen(s3_prefix) + 36;  // Estimate buffer size
    char s3_fname[s3_fname_buf_size];
    int written = snprintf(s3_fname, s3_fname_buf_size,
                           "%s_chunk%06d.fastq.zst", s3_prefix, chunk_index);

    // Check for snprintf errors / truncation
    if (written < 0 || written >= s3_fname_buf_size) {
      fprintf(stderr,
              "Child %d: Error formatting S3 path (written=%d, size=%d). "
              "Prefix or index too long?\n",
              getpid(), written, s3_fname_buf_size);
      exit(1);  // Exit child on formatting error
    }

    fprintf(stderr, "Child %d: Launching script for %s -> %s\n", getpid(),
            local_fname, s3_fname);
    execlp(script_path, script_path, local_fname, s3_fname, zstd_threads_str,
           (char *)NULL);

    // If execlp returns, it failed
    fprintf(stderr, "Child %d: Failed to exec script '%s': %s\n", getpid(),
            script_path, strerror(errno));
    exit(1);
  } else {       // Parent process
    return pid;  // Return child PID
  }
}

/**
 * @brief Opens a FASTQ file, handling .gz compression via popen/zcat. Uses
 * larger buffer for command.
 *
 * @param path Path to the FASTQ file (.gz or plain).
 * @param is_pipe Output pointer; set to 1 if popen was used (file is gzipped),
 * 0 otherwise.
 * @return FILE* A file pointer (either from fopen or popen) ready for reading.
 * Exits fatally on failure.
 */
FILE *open_fastq_file(const char *path, int *is_pipe) {
  FILE *f;
  if (strstr(path, ".gz")) {
    // Use larger buffer for command string instead of malloc, assuming 2k +
    // path is sufficient
    char cmd[2048 + strlen(path)];
    int written = snprintf(cmd, sizeof(cmd), "zcat %s", path);
    if (written < 0 || (size_t)written >= sizeof(cmd)) {
      fprintf(stderr,
              "Error: Failed to format zcat command (path too long?): %s\n",
              path);
      exit(1);  // Exit if command formatting fails
    }
    f = popen(cmd, "r");
    if (!f) {
      fprintf(stderr, "Error: popen failed for command '%s': %s\n", cmd,
              strerror(errno));
      exit(1);
    }
    *is_pipe = 1;
  } else {
    f = fopen(path, "r");
    if (!f) {
      fprintf(stderr, "Error opening input file %s: %s\n", path,
              strerror(errno));
      exit(1);
    }
    *is_pipe = 0;
  }
  return f;
}

/**
 * @brief Closes a file pointer previously opened by open_fastq_file.
 * Handles pipe streams (pclose) vs regular files (fclose) and checks return
 * status.
 *
 * @param f The file pointer to close. Can be NULL (no-op).
 * @param is_pipe Flag indicating if the file pointer originated from popen (1)
 * or fopen (0).
 */
void close_fastq_file(FILE *f, int is_pipe) {
  if (!f) return;
  if (is_pipe) {
    int status = pclose(f);
    if (status == -1) {
      perror("pclose failed");
    }  // Log pclose errors
    else if (WIFEXITED(status) && WEXITSTATUS(status) != 0) {
      fprintf(stderr,
              "Warning: Input pipe command (zcat?) exited with status %d\n",
              WEXITSTATUS(status));
    } else if (WIFSIGNALED(status)) {
      fprintf(stderr,
              "Warning: Input pipe command (zcat?) killed by signal %d\n",
              WTERMSIG(status));
    }
  } else {
    if (fclose(f) != 0) {
      perror("fclose failed");
    }  // Log fclose errors
  }
}

/**
 * @brief Reads one complete FASTQ record (4 lines) from the file stream.
 * Uses the provided buffer struct, reallocating buffers via getline as needed.
 * Performs basic FASTQ format validation. Exits fatally on error or corruption.
 * Does not support FASTQ with hard wrapping.
 *
 * @param f Input file stream.
 * @param buffers Pointer to FastqRecordBuffers struct where line pointers and
 * sizes are stored/updated.
 * @param read_num Index of record being read; only used for debug messages.
 * @return int 1 on successful read of a full record, 0 on reaching EOF cleanly.
 */
int read_fastq_record(FILE *f, FastqRecordBuffers *buffers, long read_num) {
  ssize_t title_len = getline(&(buffers->title), &(buffers->title_size), f);
  ssize_t seq_len = getline(&(buffers->seq), &(buffers->seq_size), f);
  ssize_t plus_len = getline(&(buffers->plus), &(buffers->plus_size), f);
  ssize_t quality_len =
      getline(&(buffers->quality), &(buffers->quality_size), f);

  if (title_len == -1 || seq_len == -1 || plus_len == -1 || quality_len == -1) {
    if (feof(f)) {
      if (title_len != -1 || seq_len != -1 || plus_len != -1 ||
          quality_len != -1) {
        fprintf(stderr,
                "Warning: Incomplete FASTQ record at end of file (Read pair "
                "#%ld).\n",
                read_num);
      }
      return 0;
    } else {
      perror("getline error reading FASTQ record");
      fprintf(stderr, "Error occurred reading data for read pair #%ld.\n",
              read_num);
      exit(1);
    }
  }

  if ((buffers->title)[0] != '@' || (buffers->plus)[0] != '+' ||
      seq_len != quality_len || title_len > MAX_READ_LEN ||
      seq_len > MAX_READ_LEN || plus_len > MAX_READ_LEN ||
      quality_len > MAX_READ_LEN) {
    fprintf(stderr, "Corrupt fastq detected processing read pair #%ld\n",
            read_num);
    fprintf(stderr,
            "Record starts: Title='%c' Plus='%c'. Lengths: seq=%zd qual=%zd\n",
            (buffers->title)[0], (buffers->plus)[0], seq_len, quality_len);
    exit(1);
  }
  return 1;
}
