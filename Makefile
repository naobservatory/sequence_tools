CC = gcc
CFLAGS = -Wall -Wextra -O3
LDFLAGS =

LIBS_SIZ2FASTQ = -lzstd -lz
SRC_SIZ2FASTQ = siz2fastq.c
TARGET_SIZ2FASTQ = siz2fastq

SRC_SPLIT_INTERLEAVE = split_interleave_fastqs.c
TARGET_SPLIT_INTERLEAVE = split_interleave_fastqs

TARGETS = $(TARGET_SIZ2FASTQ) $(TARGET_SPLIT_INTERLEAVE)

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin) # macOS-specific settings
    LDFLAGS += -L/usr/local/lib -L/opt/homebrew/lib
    CFLAGS += -I/usr/local/include -I/opt/homebrew/include
    BREW_PREFIX := $(shell /opt/homebrew/bin/brew --prefix)
    LDFLAGS += -L$(BREW_PREFIX)/lib
    CFLAGS += -I$(BREW_PREFIX)/include
endif

all: $(TARGETS)

$(TARGET_SIZ2FASTQ): $(SRC_SIZ2FASTQ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS_SIZ2FASTQ)

$(TARGET_SPLIT_INTERLEAVE): $(SRC_SPLIT_INTERLEAVE)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(TARGETS)

install: $(TARGETS)
	install -d $(DESTDIR)/usr/local/bin/
	install -m 755 $(TARGETS) $(DESTDIR)/usr/local/bin/

uninstall:
	rm -f $(addprefix $(DESTDIR)/usr/local/bin/,$(TARGETS))

.PHONY: all clean install uninstall

