vpath %.c ..
vpath %.c lib
vpath %.h ..
vpath %.h lib

P= runtests

OBJECTS= libunittest.so unittests_bwt.o divsufsort.o
SOURCES= bwt.c

CC= gcc
INCLUDES= -I.. -Ilib
COVERAGE= -fprofile-arcs -ftest-coverage
PROFILE= -pg
CFLAGS= -std=gnu99 -g -Wall -O0 $(INCLUDES) $(COVERAGE) $(PROFILE)
LDLIBS= -L. -Wl,-rpath,. -lunittest -lz -lm -lpthread
# Use different flags on Linux and MacOS.
ifeq ($(shell uname -s),Darwin)
	libflag= -dynamiclib
else
	libflag= -shared
endif

$(P): $(OBJECTS) $(SOURCES) runtests.c
	$(CC) $(CFLAGS) runtests.c $(OBJECTS) $(LDLIBS) -o $@

clean:
	rm -rf $(P) $(OBJECTS) \
		*.gcda *.gcno *.gcov gmon.out \
		.inspect.gdb libunittest.so.dSYM runtests.dSYM

libunittest.so: unittest.c
	$(CC) -fPIC $(libflag) $(CFLAGS) -o libunittest.so lib/unittest.c

test: $(P)
	./$(P)

inspect: $(P)
	gdb --command=.inspect.gdb --args $(P)

debug: $(P)
	gdb --command=.debug.gdb --args $(P)

valgrind: $(P)
	valgrind --leak-check=full --show-leak-kinds=all ./$(P)

vgdb: $(P)
	valgrind --vgdb=yes --vgdb-error=0 ./$(P)

unittests_mem_seed_prob.gcda gmon.out: $(P)
	-./$(P)

cov: unittests_mem_seed_prob.gcda
	gcov unittests_mem_seed_prob.c
	cat mem_seed_prob.c.gcov

prof: gmon.out
	gprof ./$(P) --flat-profile --brief
