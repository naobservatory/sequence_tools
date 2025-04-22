#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zstd.h>
#include <zlib.h>
#include <unistd.h>

#define BUFFER_SIZE 4096
#define LINE_BUFFER_SIZE 1024 * 1024

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
        output.regular_file = NULL;
    } else {
        output.regular_file = fopen(filename, "wb");
        output.gz_file = NULL;
    }

    return output;
}

size_t write_output(OutputFile* file, const void* buffer, size_t size) {
    if (file->is_gzipped) {
        return gzwrite(file->gz_file, buffer, size);
    } else {
        return fwrite(buffer, /*nitems=*/1, size, file->regular_file);
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

    char in_buffer[BUFFER_SIZE];

    ZSTD_DStream* dstream = ZSTD_createDStream();
    if (dstream == NULL) {
        fprintf(stderr, "Error creating zstd decompression stream\n");
        return 1;
    }

    ZSTD_inBuffer input = { in_buffer, 0, 0 };

    char out_buffer[BUFFER_SIZE];
    char line_buffer[LINE_BUFFER_SIZE];
    size_t line_pos = 0;

    OutputFile out_file1 = open_output_file(argv[optind], use_gzip);
    OutputFile out_file2 = open_output_file(argv[optind + 1], use_gzip);

    // Verify files opened successfully.
    if ((!out_file1.gz_file && !out_file1.regular_file) ||
        (!out_file2.gz_file && !out_file2.regular_file)) {
        fprintf(stderr, "Error opening output files\n");
        ZSTD_freeDStream(dstream);
        return 1;
    }

    long long line_count = 0;
    OutputFile* current_file = &out_file1;

    size_t read_size;
    while ((read_size = fread(
              in_buffer, /*nitems=*/1, BUFFER_SIZE, stdin)) > 0) {
        input.size = read_size;
        input.pos = 0;

        while (input.pos < input.size) {
            ZSTD_outBuffer output = { out_buffer, BUFFER_SIZE, 0 };

            size_t result = ZSTD_decompressStream(dstream, &output, &input);
            if (ZSTD_isError(result)) {
                fprintf(stderr,
                        "Error during decompression: %s\n",
                        ZSTD_getErrorName(result));
                close_output_file(&out_file1);
                close_output_file(&out_file2);
                ZSTD_freeDStream(dstream);
                return 1;
            }

            for (size_t i = 0; i < output.pos; i++) {
                char c = out_buffer[i];
                line_buffer[line_pos++] = c;

                if (c == '\n') {
                    write_output(current_file, line_buffer, line_pos);
                    line_pos = 0;

                    line_count++;
                    if (line_count % 8 < 4) {
                        current_file = &out_file1;
                    } else {
                        current_file = &out_file2;
                    }
                }

                if (line_pos >= LINE_BUFFER_SIZE - 1) {
                  write_output(current_file, line_buffer, line_pos);
                    line_pos = 0;
                }
            }
        }
    }

    // Maybe we have some left over, if the file doesn't end with newline.
    if (line_pos > 0) {
        write_output(current_file, line_buffer, line_pos);
    }

    close_output_file(&out_file1);
    close_output_file(&out_file2);
    ZSTD_freeDStream(dstream);

    return 0;
}
