# Autoconf/automake file

bin_PROGRAMS = test_cvmhlabn unittest

# General compiler/linker flags
AM_CFLAGS = -Wall -O3 -std=c99 -D_LARGEFILE_SOURCE \
        -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 ${CFLAGS} -I../src
#AM_LDFLAGS = ${LDFLAGS} -L../src -lcvmhlabn
AM_LDFLAGS = ${LDFLAGS} -L../src -lcvmhlabn -L../gctpc/source -lgeo -lm

# Dist sources
unittest_SOURCES = *.c *.h

TARGETS = $(bin_PROGRAMS)

.PHONY = run_api run_unit

all: $(bin_PROGRAMS)

unittest: unittest.o unittest_defs.o test_helper.o test_vx_sub.o \
	test_vx_exec.o test_vx_lite_exec.o
	$(CC) -o $@ $^ $(AM_LDFLAGS)

test_cvmhlabn: test_cvmhlabn.o
	$(CC) -g -o $@ $^ $(AM_CFLAGS) $(AM_LDFLAGS) 

run_api: test_cvmhlabn
	./run_test_cvmhlabn.sh

run_unit : unittest
	./unittest

clean:
	rm -rf *~ *.o $(bin_PROGRAMS) 

install:
	mkdir -p ${prefix}/tests
