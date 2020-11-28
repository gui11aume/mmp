P= mmp

CC= gcc
CFLAGS= -std=c99 -Wall
OBJECTS= bwt.o map.o QSufSort.o bwt_gen.o sesame.o

# development flags: -DNOQUAL -DFASTOUT
all: CFLAGS += -DNDEBUG -O3
all: $(P)

debug: CFLAGS += -DDEBUG -g -O0
debug: $(P)

profile: CFLAGS += -pg -O0 -DNOQUAL
profile: $(P)

analyze: CC= clang --analyze
analyze: CFLAGS += -DDEBUG -g -O0
analyze: $(OBJECTS)


$(P): mmp.c bwt.h $(OBJECTS)
	$(CC) $(CFLAGS) mmp.c $(OBJECTS) -o $(P) -lm -lpthread
#	$(CC) $(CFLAGS) mmp.c $(OBJECTS) -o $(P) -Wl,--no-as-needed -lprofiler -lm

clean:
	rm -f $(OBJECTS) $(P)
