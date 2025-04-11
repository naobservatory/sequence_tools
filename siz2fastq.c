#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zstd.h>
#include <zlib.h>

#define BUFFER_SIZE 4096
#define LINE_BUFFER_SIZE 1024 * 1024

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr,
            "Usage: %s output_1.fastq.gz output_2.fastq.gz\n",
            argv[0]);
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
  
  gzFile out_file1 = gzopen(argv[1], "wb");
  gzFile out_file2 = gzopen(argv[2], "wb");
  if (!out_file1 || !out_file2) {
    fprintf(stderr, "Error opening output files\n");
    ZSTD_freeDStream(dstream);
    return 1;
  }
    
  int line_count = 0;
  gzFile current_file = out_file1;
    
  size_t read_size;
  while ((read_size = fread(in_buffer, 1, BUFFER_SIZE, stdin)) > 0) {
    input.size = read_size;
    input.pos = 0;
        
    while (input.pos < input.size) {
      ZSTD_outBuffer output = { out_buffer, BUFFER_SIZE, 0 };
            
      size_t result = ZSTD_decompressStream(dstream, &output, &input);
      if (ZSTD_isError(result)) {
        fprintf(stderr,
                "Error during decompression: %s\n",
                ZSTD_getErrorName(result));
        gzclose(out_file1);
        gzclose(out_file2);
        ZSTD_freeDStream(dstream);
        return 1;
      }
            
      for (size_t i = 0; i < output.pos; i++) {
        char c = out_buffer[i];

        line_buffer[line_pos++] = c;
                
        if (c == '\n') {
          gzwrite(current_file, line_buffer, line_pos);
          line_pos = 0;
                    
          line_count++;
          if (line_count % 8 < 4) {
            current_file = out_file1;
          } else {
            current_file = out_file2;
          }
        }
                
        if (line_pos >= LINE_BUFFER_SIZE - 1) {
          gzwrite(current_file, line_buffer, line_pos);
          line_pos = 0;
        }
      }
    }
  }

  // Maybe we have some left over, if the file doesn't end with newline.
  if (line_pos > 0) {
    gzwrite(current_file, line_buffer, line_pos);
  }
    
  gzclose(out_file1);
  gzclose(out_file2);
  ZSTD_freeDStream(dstream);
  
  return 0;
}
