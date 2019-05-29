P= smmfdp

CC= gcc
CFLAGS= -std=c99 -Wall
OBJECTS= bwt.o map.o divsufsort.o sesame.o

# development flags: -DNOQUAL -DFASTOUT
all: CFLAGS += -DNDEBUG -O3
all: $(P)

debug: CFLAGS += -DDEBUG -g -O0
debug: $(P)

profile: CFLAGS += -pg -O0 -DNOQUAL
profile: $(P)

noqual: CFLAGS += -DNDEBUG -DNOQUAL -O3
noqual: $(P)

$(P): smmfdp.c bwt.h $(OBJECTS)
	$(CC) $(CFLAGS) smmfdp.c $(OBJECTS) -o $(P) -lm

clean:
	rm -f $(OBJECTS) $(P)
