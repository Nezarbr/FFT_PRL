CC=gcc
CFLAGS=-fopenmp -Wall -O3 -mavx2 -lm

SRCS := $(wildcard *.c)
BINS := $(patsubst %.c,%,$(SRCS))

all: $(BINS)

%: %.c
	$(CC) -o $@ $< $(CFLAGS) 

clean:
	rm -f $(BINS)