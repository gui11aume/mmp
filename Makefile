P= smmfdp

CC= gcc
CFLAGS= -std=c99 -Wall
OBJECTS= bwt.o divsufsort.o

all: CFLAGS += -DNDEBUG -O3
all: $(P)

debug: CFLAGS += -DDEBUG -g -O0
debug: $(P)

$(P): smmfdp.c divsufsort.o bwt.o bwt.h
	$(CC) $(CFLAGS) smmfdp.c $(OBJECTS) -o $(P)

clean:
	rm -f $(OBJECTS) $(P)
