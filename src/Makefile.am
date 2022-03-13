# Autoconf/automake file

lib_LIBRARIES = libvxapi.a libvxapi_cvmhlabn.a libcvmhlabn.a 
bin_PROGRAMS = vx vx_lite vx_slice vx_lite_cvmhlabn vx_cvmhlabn run_vx.sh run_vx_lite.sh
include_HEADERS = vx_sub.h

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS} -Wall -O3 -std=c99 -D_LARGEFILE_SOURCE \
            -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -fPIC
AM_LDFLAGS = ${LDFLAGS} -L../gctpc/source -lgeo

# Dist sources
libcvmhlabn_a_SOURCES = vx_sub_cvmhlabn.c vx_io.c 
libvxapi_a_SOURCES = vx_sub.c scec1d.c vs30_gtl.c vx_io.c utils.c 
vx_lite_cvmhlabn_SOURCES = vx_lite_cvmhlabn.c
vx_cvmhlabn_SOURCES = vx_cvmhlabn.c
vx_slice_SOURCES = vx_slice.c
vx_lite_SOURCES = vx_lite.c
run_vx_sh_SOURCES = run_vx.sh
run_vx_lite_sh_SOURCES = run_vx_lite.sh

TARGETS = vx_lite_cvmhlabn vx_cvmhlabn libvxapi_cvmhlabn.a libvxapi.a libcvmhlabn.a libcvmhlabn.so

all: $(TARGETS)

install:
	mkdir -p ${prefix}
	mkdir -p ${prefix}/lib
	mkdir -p ${prefix}/include
	cp libcvmhlabn.so ${prefix}/lib
	cp libcvmhlabn.a ${prefix}/lib
	cp cvmhlabn.h ${prefix}/include

libcvmhlabn.a: vx_sub_cvmhlabn.o vx_io.o utils.o cvmhlabn_static.o 
	$(AR) rcs $@ $^

cvmhlabn_static.o: cvmhlabn.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)

libcvmhlabn.so: vx_sub_cvmhlabn.o vx_io.o utils.o cvmhlabn.o
	$(CC) -shared $(AM_CFLAGS) -o libcvmhlabn.so $^ $(AM_LDFLAGS)

libvxapi.a: vx_sub.o scec1d.o vs30_gtl.o vx_io.o utils.o *.h
	$(AR) rcs $@ $^

libvxapi_cvmhlabn.a: vx_sub_cvmhlabn.o vx_io.o utils.o *.h
	$(AR) rcs $@ $^

cvmhlabn.o: cvmhlabn.c
	$(CC) -fPIC -DDYNAMIC_LIBRARY -o $@ -c $^ $(AM_CFLAGS)

vx_lite_cvmhlabn.o : vx_lite_cvmhlabn.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)

vx_lite_cvmhlabn : vx_lite_cvmhlabn.o libvxapi_cvmhlabn.a
	$(CC) -o $@ $^ $(AM_LDFLAGS)

vx_cvmhlabn.o : vx_cvmhlabn.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)

vx_cvmhlabn : vx_cvmhlabn.o libcvmhlabn.a
	$(CC) -o $@ $^ $(AM_LDFLAGS)

vx: vx.o
	$(CC) -o $@ $^ $(AM_LDFLAGS)

vx_lite: vx_lite.o libvxapi.a
	$(CC) -o $@ $^ $(AM_LDFLAGS)

vx_slice: vx_slice.o libvxapi.a
	$(CC) -o $@ $^ $(AM_LDFLAGS)

run_vx.sh:
run_vx_lite.sh:


clean:
	rm -rf $(TARGETS)
	rm -rf *.o 

