CC = gcc
CFLAGS = -Wall -Wextra -O2
LDFLAGS =
LIBS = -lzstd -lz
TARGET = siz2fastq
SRC = siz2fastq.c

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin) # macOS-specific settings
    LDFLAGS += -L/usr/local/lib -L/opt/homebrew/lib
    CFLAGS += -I/usr/local/include -I/opt/homebrew/include
    BREW_PREFIX := $(shell /opt/homebrew/bin/brew --prefix)
    LDFLAGS += -L$(BREW_PREFIX)/lib
    CFLAGS += -I$(BREW_PREFIX)/include
endif

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

clean:
	rm -f $(TARGET)

install: $(TARGET)
	install -d $(DESTDIR)/usr/local/bin/
	install -m 755 $(TARGET) $(DESTDIR)/usr/local/bin/

uninstall:
	rm -f $(DESTDIR)/usr/local/bin/$(TARGET)

.PHONY: all clean install uninstall
