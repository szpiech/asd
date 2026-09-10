
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
# No -Wall for the vendored netlib sources; they are not ours to fix and the noise
#  hides warnings from src/.
FFLAGS=-c $(OPT) -g
LFLAGS=-g -o

bin/lodestar: src/Main.o
	mkdir -p bin
	$(CC) $(LFLAGS) bin/lodestar src/*.o lib/kstring.o lib/lapack/*.o lib/gsl/*.o -lz -lm -lpthread -lgfortran

# NOTE: each object now lists its own source and headers as prerequisites. Previously the
#  rules named only other objects, so editing a .c or .h file did not trigger a rebuild
#  and `make` could quietly link stale objects.
src/Main.o: src/Main.c src/LODESTAR.h src/Interface.h src/LODESTAR.o src/Interface.o
	$(CC) $(CFLAGS) src/Main.c -o src/Main.o

src/Interface.o: src/Interface.c src/Interface.h lib/kstring.o
	$(CC) $(CFLAGS) src/Interface.c -o src/Interface.o

src/LODESTAR.o: src/LODESTAR.c src/LODESTAR.h src/BlockList.h lib/gsl src/HaplotypeEncoder.o src/BlockList.o src/MatrixOperations.o
	$(CC) $(CFLAGS) -DHAVE_INLINE src/LODESTAR.c -o src/LODESTAR.o

src/BlockList.o: src/BlockList.c src/BlockList.h
	$(CC) $(CFLAGS) src/BlockList.c -o src/BlockList.o

src/HaplotypeEncoder.o: src/HaplotypeEncoder.c src/HaplotypeEncoder.h src/VCFLocusParser.o
	$(CC) $(CFLAGS) src/HaplotypeEncoder.c -o src/HaplotypeEncoder.o

src/VCFLocusParser.o: src/VCFLocusParser.c src/VCFLocusParser.h
	$(CC) $(CFLAGS) src/VCFLocusParser.c -o src/VCFLocusParser.o

src/MatrixOperations.o: src/MatrixOperations.c src/MatrixOperations.h lib/lapack 
	$(CC) $(CFLAGS) src/MatrixOperations.c -o src/MatrixOperations.o

# NOTE: this used to write its output to src/kstring.o, so the target file never existed
#  and the rule was re-run on every invocation of make.
lib/kstring.o:
	$(CC) $(CFLAGS) lib/kstring.c -o lib/kstring.o

# NOTE: the old rule compiled into the repository root and then `mv *.o lib/lapack`,
#  which races with the C rules under `make -j` and can move their objects too. Compiling
#  each source to an explicit -o keeps the rule parallel-safe. It also used $(CFLAGS),
#  which passed the C-only `-I lib` to gfortran.
.PHONY: lib/lapack
lib/lapack:
	@for f in lib/lapack/*.f; do \
		$(FC) $(FFLAGS) $$f -o $${f%.f}.o || exit 1; \
	done

.PHONY: lib/gsl
lib/gsl:
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/error.c -o lib/gsl/error.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/message.c -o lib/gsl/message.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/stream.c -o lib/gsl/stream.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/default.c -o lib/gsl/default.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/rng.c -o lib/gsl/rng.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/mt.c -o lib/gsl/mt.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/types.c -o lib/gsl/types.o


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