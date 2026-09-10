
# File: Makefile
# Date: 9 May 2024
# Author: T. Quinn Smith
# Principal Investigator: Dr. Zachary A. Szpiech
# Purpose: Build LODESTAR. 

# NOTE: these are `=` and not `?=`. GNU make defines CC = cc and FC = f77 as built-in
#  defaults, and `?=` only assigns to a variable that is UNDEFINED -- a built-in default
#  counts as defined. With `?=` a bare `make` therefore invoked f77 (which does not exist
#  on most systems today, so the build failed at the LAPACK step) and compiled the C with
#  cc, which on macOS is clang, contradicting the gfortran requirement in
#  MatrixOperations.c. Override on the command line as usual: make CC=clang.
CC=gcc
FC=gfortran

# -O2 matters a great deal here: at the previous flags (-g with no -O) every C source and
#  all 99 vendored LAPACK/BLAS Fortran files were compiled unoptimised, which cost a
#  factor of 2.4-3.8 in wall-clock time on the same source.
OPT=-O2
CFLAGS=-c -Wall $(OPT) -g -I lib
# FOPT is kept separate from OPT so the sanitizer build can instrument the C sources
#  without instrumenting the vendored Fortran. Sanitizing 99 LAPACK files finds nothing
#  (the bugs are in src/) and requires the C and Fortran compilers to ship the same
#  sanitizer runtime version, which is not true for a mixed clang/gfortran toolchain.
FOPT=-O2
# No -Wall for the vendored netlib sources; they are not ours to fix and the noise
#  hides warnings from src/.
FFLAGS=-c $(FOPT) -g
LFLAGS=-g -o

# LIBDIRS is a hook for -L paths. It is needed whenever the C compiler does not already
#  know where libgfortran lives -- notably clang on macOS with a Homebrew gfortran:
#    make CC=clang LIBDIRS="-L$$(brew --prefix gcc)/lib/gcc/current"
#  With CC=gcc and a matching gfortran, the compiler finds it and LIBDIRS can stay empty.
LIBDIRS=
LIBS=-lz -lm -lpthread -lgfortran

# Object lists, so every prerequisite below is a real file. The vendored libraries used to
#  be built by .PHONY targets named lib/lapack and lib/gsl, which meant they ran on every
#  invocation and -- because src/LODESTAR.o and src/MatrixOperations.o depended on them --
#  forced a full recompile and relink of the project on an unchanged tree.
SRC_OBJS=src/VCFLocusParser.o src/HaplotypeEncoder.o src/BlockList.o \
         src/MatrixOperations.o src/LODESTAR.o src/Interface.o src/Main.o
GSL_OBJS=lib/gsl/error.o lib/gsl/message.o lib/gsl/stream.o lib/gsl/default.o \
         lib/gsl/rng.o lib/gsl/mt.o lib/gsl/types.o
LAPACK_OBJS=$(patsubst %.f,%.o,$(wildcard lib/lapack/*.f))
OBJS=$(SRC_OBJS) lib/kstring.o $(LAPACK_OBJS) $(GSL_OBJS)

bin/lodestar: $(OBJS)
	mkdir -p bin
	$(CC) $(LFLAGS) bin/lodestar $(OBJS) $(LIBDIRS) $(LIBS)

# NOTE: each object now lists its own source and headers as prerequisites. Previously the
#  rules named only other objects, so editing a .c or .h file did not trigger a rebuild
#  and `make` could quietly link stale objects.
src/Main.o: src/Main.c src/LODESTAR.h src/Interface.h src/MatrixOperations.h
	$(CC) $(CFLAGS) src/Main.c -o src/Main.o

src/Interface.o: src/Interface.c src/Interface.h
	$(CC) $(CFLAGS) src/Interface.c -o src/Interface.o

src/LODESTAR.o: src/LODESTAR.c src/LODESTAR.h src/BlockList.h src/VCFLocusParser.h src/HaplotypeEncoder.h src/MatrixOperations.h
	$(CC) $(CFLAGS) -DHAVE_INLINE src/LODESTAR.c -o src/LODESTAR.o

src/BlockList.o: src/BlockList.c src/BlockList.h
	$(CC) $(CFLAGS) src/BlockList.c -o src/BlockList.o

src/HaplotypeEncoder.o: src/HaplotypeEncoder.c src/HaplotypeEncoder.h src/VCFLocusParser.h
	$(CC) $(CFLAGS) src/HaplotypeEncoder.c -o src/HaplotypeEncoder.o

src/VCFLocusParser.o: src/VCFLocusParser.c src/VCFLocusParser.h
	$(CC) $(CFLAGS) src/VCFLocusParser.c -o src/VCFLocusParser.o

src/MatrixOperations.o: src/MatrixOperations.c src/MatrixOperations.h
	$(CC) $(CFLAGS) src/MatrixOperations.c -o src/MatrixOperations.o

# NOTE: this used to write its output to src/kstring.o, so the target file never existed,
#  the rule re-ran on every invocation of make, and the object was only picked up by the
#  link because `src/*.o` happened to glob it.
lib/kstring.o: lib/kstring.c lib/kstring.h
	$(CC) $(CFLAGS) lib/kstring.c -o lib/kstring.o

# NOTE: the LAPACK rule used to compile into the repository root and then `mv *.o
#  lib/lapack`, which races with the C rules under `make -j` and can move their objects
#  too. It also passed $(CFLAGS), handing the C-only `-I lib` to gfortran. A pattern rule
#  removes both problems and makes the build correctly incremental and -j safe.
lib/lapack/%.o: lib/lapack/%.f
	$(FC) $(FFLAGS) $< -o $@

lib/gsl/%.o: lib/gsl/%.c
	$(CC) $(CFLAGS) -DHAVE_INLINE $< -o $@


# A build with the address and undefined-behaviour sanitizers, for use before a release.
#  Run `make asan && ./bin/lodestar_asan ...`.
.PHONY: asan
asan:
	$(MAKE) clean
	$(MAKE) OPT="-O1 -fsanitize=address,undefined -fno-omit-frame-pointer" \
	        LFLAGS="-g -fsanitize=address,undefined -o"
	mv bin/lodestar bin/lodestar_asan

.PHONY: test
test: bin/lodestar
	./tests/run_tests.sh

PREFIX?=/usr/local
.PHONY: install
install: bin/lodestar
	mkdir -p $(PREFIX)/bin
	cp bin/lodestar $(PREFIX)/bin/lodestar

# NOTE: -f, so `make clean` on an already-clean tree succeeds instead of exiting 1.
.PHONY: clean
clean:
	rm -f bin/lodestar bin/lodestar_asan src/*.o lib/lapack/*.o lib/*.o lib/gsl/*.o