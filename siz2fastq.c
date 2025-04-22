#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zstd.h>
#include <zlib.h>
#include <unistd.h>

#define BUFFER_SIZE 4096
#define LINE_BUFFER_SIZE 1024 * 1024

// Split an interleaved zstd-compressed fastq file into a pair of optionally
// gzip-compressed fastq files.  No input validation is performed: the input
// must be valid fastq with no hard wrapping or you'll get garbage output.

// Struct to track what kind of file we're writing, since we optionally
// compress our output.
typedef struct {
    FILE* regular_file;
    gzFile gz_file;
    bool is_gzipped;
} OutputFile;

OutputFile open_output_file(const char* filename, bool use_gzip) {
    OutputFile output;
    output.is_gzipped = use_gzip;

    if (use_gzip) {
        output.gz_file = gzopen(filename, "wb");
        if (output.gz_file == Z_NULL) {
            fprintf(stderr,
                    "Unable to open %s as a gzip-compressed output file\n",
                    filename);
            exit(1);
        }
        output.regular_file = NULL;
    } else {
        output.regular_file = fopen(filename, "wb");
        if (output.regular_file == NULL) {
            fprintf(stderr, "Unable to open %s as output file\n", filename);
            exit(1);
        }
        output.gz_file = NULL;
    }

    return output;
}

// Write `size` bytes from `buffer` to `file` or exit with an error message.
void write_output_or_die(
         OutputFile* file, const void* buffer, size_t size) {
    size_t bytes_written;
    if (file->is_gzipped) {
        bytes_written = gzwrite(file->gz_file, buffer, size);
    } else {
        bytes_written = fwrite(buffer, /*nitems=*/1, size, file->regular_file);
    }
    if (bytes_written != size) {
        fprintf(stderr, "Failed to write all data\n");
        exit(1);
    }
}

void close_output_file(OutputFile* file) {
    if (file->is_gzipped) {
        gzclose(file->gz_file);
    } else {
        fclose(file->regular_file);
    }
}

void print_usage(const char* program_name) {
    fprintf(stderr,
            "Usage: %s [-z] out1.fastq[.gz] out2.fastq[.gz]\n",
            program_name);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -z    Compress output with gzip\n");
}

int main(int argc, char *argv[]) {
    bool use_gzip = false;
    int opt;

    // Parse command line options.
    while ((opt = getopt(argc, argv, "z")) != -1) {
        switch (opt) {
            case 'z':
                use_gzip = true;
                break;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    // Check if we have the correct number of arguments after options.  The way
    // getopt works is a bit weird, and it uses some global variables (optarg,
    // opterr, optind) for input and output.  Here we're using optind, which is
    // the index of the next unprocessed argument.  Verify we have exactly two
    // arguments left (one for each output file) after getopt processes our
    // flags.
    if (argc - optind != 2) {
        print_usage(argv[0]);
        return 1;
    }

    // Unprocessed bytes we've read from the input.
    char out_buffer[BUFFER_SIZE];
    // Unprocessed bytes of the line we're currently reading.
    char line_buffer[LINE_BUFFER_SIZE];
    // How far we are into the line we're currently reading.
    size_t line_pos = 0;

    // As discussed above, optind is the index of the first positional
    // argument.
    OutputFile out1 = open_output_file(argv[optind], use_gzip);
    OutputFile out2 = open_output_file(argv[optind + 1], use_gzip);

    char in_buffer[BUFFER_SIZE];

    ZSTD_DStream* dstream = ZSTD_createDStream();
    if (dstream == NULL) {
        fprintf(stderr, "Error creating zstd decompression stream\n");
        return 1;
    }

    ZSTD_inBuffer input = { in_buffer, 0, 0 };

    // The number of lines we've read so far.
    long long line_count = 0;
    // We switch back and forth between writing to out1 and out2, four lines to
    // one then four lines to the other.
    OutputFile* current_file = &out1;

    // Loop reading from stdin until there's no more input.  We read bytes not
    // lines, so we do need to handle cases where a read gives us multiple
    // lines or a partial line.
    size_t read_size;
    while ((read_size = fread(
              in_buffer, /*nitems=*/1, BUFFER_SIZE, stdin)) > 0) {
        input.size = read_size;
        input.pos = 0;

        // In most cases the call to ZSTD_decompressStream below will consume
        // all of `input`, but not it:
        //  1. `input` represents multiple concatenated zstd files, in which
        //     case it will stop partway through the input.  When it does that
        //     it sets input.pos to how far it got, and we'll just call again
        //     on the remainder because we don't care that it's concatenated.
        //  2. We can't fit all of the decompressed data in the output buffer,
        //     in which case we'll consume what we have and loop to handle the
        //     rest.
        while (input.pos < input.size) {
            ZSTD_outBuffer output = { out_buffer, BUFFER_SIZE, 0 };

            size_t result = ZSTD_decompressStream(dstream, &output, &input);
            if (ZSTD_isError(result)) {
                fprintf(stderr,
                        "Error during decompression: %s\n",
                        ZSTD_getErrorName(result));
                close_output_file(&out1);
                close_output_file(&out2);
                ZSTD_freeDStream(dstream);
                return 1;
            }

            for (size_t i = 0; i < output.pos; i++) {
                // Save the input characters one at a time to line_buffer,
                // building up the current line.
                char c = out_buffer[i];
                line_buffer[line_pos++] = c;

                // When we learn the line has ended, process and reset the
                // line buffer.
                if (c == '\n') {
                    write_output_or_die(current_file, line_buffer, line_pos);

                    // Reset the buffer, by telling ourselves it's empty.
                    line_pos = 0;

                    // Every four lines (title, seq, plus line, quality) we
                    // switch which file to write to, deinterleaving.
                    line_count++;
                    if (line_count % 8 < 4) {
                        current_file = &out1;
                    } else {
                        current_file = &out2;
                    }
                }

                // If we've built up a very long line in process, flush it out.
                // Since we don't need to know how far we are into the line
                // this isn't a problem (it would be if we were doing more
                // involved processing that deinterleaving).  Note that this is
                // very unlikely to happen, unless we have interleaved
                // long-read data,
                if (line_pos >= LINE_BUFFER_SIZE - 1) {
                    write_output_or_die(current_file, line_buffer, line_pos);
                    line_pos = 0;
                }
            }
        }
    }

    // Maybe we have some left over, if the file doesn't end with newline.
    if (line_pos > 0) {
        write_output_or_die(current_file, line_buffer, line_pos);
    }

    // These would be closed automatically on exit, but it's tidier to clean
    // them up.
    close_output_file(&out1);
    close_output_file(&out2);
    ZSTD_freeDStream(dstream);

    return 0;
}
