# Autoconf/automake file

bin_PROGRAMS = cvmhlabn_api_validate cvmhlabn_vxlite_validate

# General compiler/linker flags
AM_CFLAGS = -DDYNAMIC_LIBRARY -Wall -O3 -std=c99 -D_LARGEFILE_SOURCE \
        -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -I../src
AM_LDFLAGS = -L../src -lcvmhlabn -L../gctpc/source -lgctpc -lm -ldl

# Dist sources
cvmhlabn_api_validate_SOURCES = cvmhlabn_api_validate.c
cvmhlabn_vxlite_validate_SOURCES = cvmhlabn_vxlite_validate.c

.PHONY = run_validate run_validate_with_ucvm

all: $(bin_PROGRAMS)

############################################
# Executables
############################################

cvmhlabn_vxlite_validate : cvmhlabn_vxlite_validate.o 
	$(CC) -o $@ $^ $(AM_LDFLAGS)

cvmhlabn_api_validate : cvmhlabn_api_validate.o 
	$(CC) -o $@ $^ $(AM_LDFLAGS)

### only run after ucvm is built and installed
if UCVM
cvmhlabn_ucvm_validate.o : cvmhlabn_ucvm_validate.c
	$(CC) -c -o $@ $^ $(AM_CFLAGS) $(CFLAGS)

cvmhlabn_ucvm_validate : cvmhlabn_ucvm_validate.o 
	$(CC) -o $@ $^ $(AM_LDFLAGS) $(LDFLAGS)
endif

############################################
# Special build targets
############################################

run_validate: cvmhlabn_api_validate cvmhlabn_vxlite_validate
	./run_vxlite_validate.sh
	./run_api_validate.sh

run_validate_with_ucvm: cvmhlabn_ucvm_validate
	./run_ucvm_validate.sh

clean:
	rm -rf *~ *.o $(bin_PROGRAMS) cvmhlabn_ucvm_validate *.txt

