# Autoconf/automake file

bin_PROGRAMS = test_cvmhlabn

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS} -I../src
AM_LDFLAGS = ${LDFLAGS} -L../src -lcvmhlabn

objects = test_cvmhlabn.o
TARGETS = $(bin_PROGRAMS)

.PHONY = run_api

all: $(bin_PROGRAMS)

install:
	mkdir -p ${prefix}/tests
	cp apitest ${prefix}/tests
	cp run_test_cvmhlabn.sh ${prefix}/tests

apitest: $(objects)
	$(CC) -g -o $@ $^ $(AM_CFLAGS) $(AM_LDFLAGS) 

$(objects): %.o: %.c
	$(CC) -g -o $@ -c $^ $(AM_CFLAGS) -I${prefix}/include

clean:
	rm -rf *~ *.o test_cvmlabn 

run_api: test_cvmhlabn
	./run_test_cvmhlabn.sh

