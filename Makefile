P= mmp

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

$(P): mmp.c bwt.h $(OBJECTS)
	$(CC) $(CFLAGS) mmp.c $(OBJECTS) -o $(P) -lm
#	$(CC) $(CFLAGS) mmp.c $(OBJECTS) -o $(P) -Wl,--no-as-needed -lprofiler -lm

clean:
	rm -f $(OBJECTS) $(P)
