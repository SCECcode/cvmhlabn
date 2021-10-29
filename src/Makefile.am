# Autoconf/automake file

lib_LIBRARIES = libcvmhlabn.so 
bin_PROGRAMS = vx_lite_cvmhlabn
include_HEADERS = vx_io.h vx_sub.h voxet.h coor_para.h cproj.h cvmhlabn.h params.h proj.h

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS} -Wall -O3 -std=c99 -D_LARGEFILE_SOURCE \
            -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -fPIC
AM_LDFLAGS = ${LDFLAGS} -L../gctpc/source -lgeo

# Dist sources
libcvmhlabn_so_SOURCES = vx_sub.c vx_io.c *.h
vx_lite_cvmhlabn_SOURCES = vx_lite_cvmhlabn.c

all: $(lib_LIBRARIES) $(bin_PROGRAMS)

TARGETS = vx_lite_cvmhlabn libvxapi.a libcvmhlabn.a libcvmhlabn.so

all: $(TARGETS)

install:
	mkdir -p ${prefix}
	mkdir -p ${prefix}/lib
	mkdir -p ${prefix}/include
	cp libcvmhlabn.so ${prefix}/lib
	cp libcvmhlabn.a ${prefix}/lib
	cp cvmhlabn.h ${prefix}/include

libcvmhlabn.a: vx_sub.o vx_io.o utils.o cvmhlabn_static.o 
	$(AR) rcs $@ $^

cvmhlabn_static.o: cvmhlabn.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)

libcvmhlabn.so: vx_sub.o vx_io.o utils.o cvmhlabn.o
	$(CC) -shared $(AM_CFLAGS) -o libcvmhlabn.so $^ $(AM_LDFLAGS)

libvxapi.a: vx_sub.o vx_io.o utils.o *.h
	$(AR) rcs $@ $^

cvmhlabn.o: cvmhlabn.c
	$(CC) -fPIC -DDYNAMIC_LIBRARY -o $@ -c $^ $(AM_CFLAGS)

vx_lite_cvmhlabn.o : vx_lite_cvmhlabn.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)

vx_lite_cvmhlabn : vx_lite_cvmhlabn.o libvxapi.a
	$(CC) -o $@ $^ $(AM_LDFLAGS)

run_vx.sh:

run_vx_lite.sh:


clean:
	rm -rf $(TARGETS)
	rm -rf *.o 

