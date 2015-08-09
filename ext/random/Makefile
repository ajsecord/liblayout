COMMON=-Wall -ansi -fstrict-aliasing -g

# Optimisation with assert code
CXXFLAGS=${COMMON} -O3

# Optimisation without assert code
#CXXFLAGS=${COMMON} -O3 -DNDEBUG

# Profiling with optimisations without assert code
#CXXFLAGS=${COMMON} -pg -O3 -DNDEBUG

# Debugging
#CXXFLAGS=${COMMON}

# C compiler
CC=gcc

# Copy C++ flags to C flags
CFLAGS=${CXXFLAGS}

LDFLAGS=-lm
PROGS=test_uniform test_gaussian test_exponential test_sobol test_halton \
		test_hammersley test_jittered_halton test_jittered_hammersley
TESTS=do_uniform_2d do_uniform_3d do_gaussian_2d do_gaussian_3d \
        do_exponential_2d do_exponential_3d do_sobol_2d do_sobol_3d \
		do_halton_2d do_halton_3d do_hammersley_2d do_hammersley_3d \
		do_jittered_halton_2d do_jittered_halton_3d \
		do_jittered_hammersley_2d do_jittered_hammersley_3d


SOURCES=random.c test_exponential.c test_gaussian.c test_uniform.c \
		test_sobol.c test_halton.c test_hammersley.c test_jittered_halton.c \
		test_jittered_hammersley.c

all: $(PROGS)

ifeq (.depend, $(wildcard .depend))
include .depend
endif

tests: $(PROGS) $(TESTS)

test_uniform: random.o test_uniform.c
	$(CC) $(CFLAGS) -o $@ $@.c random.o $(LDFLAGS)

test_gaussian: random.o test_gaussian.c
	$(CC) $(CFLAGS) -o $@ $@.c random.o $(LDFLAGS)

test_exponential: random.o test_exponential.c
	$(CC) $(CFLAGS) -o $@ $@.c random.o $(LDFLAGS)

test_sobol: random.o test_sobol.c
	$(CC) $(CFLAGS) -o $@ $@.c random.o $(LDFLAGS)

test_halton: random.o test_halton.c
	$(CC) $(CFLAGS) -o $@ $@.c random.o $(LDFLAGS)

test_hammersley: random.o test_hammersley.c
	$(CC) $(CFLAGS) -o $@ $@.c random.o $(LDFLAGS)

test_jittered_halton: random.o test_jittered_halton.c
	$(CC) $(CFLAGS) -o $@ $@.c random.o $(LDFLAGS)

test_jittered_hammersley: random.o test_jittered_hammersley.c
	$(CC) $(CFLAGS) -o $@ $@.c random.o $(LDFLAGS)

do_uniform_2d:
	./test_uniform 10000 2 > test.dat
	echo "plot 'test.dat' with dots; pause 5" | gnuplot

do_uniform_3d:
	./test_uniform 10002 3 > test.dat
	echo "set parametric; splot 'test.dat' with dots; pause 5" | gnuplot

do_gaussian_2d:
	./test_gaussian 10000 2 > test.dat
	echo "plot 'test.dat' with dots; pause 5" | gnuplot

do_gaussian_3d:
	./test_gaussian 10002 3 > test.dat
	echo "set parametric; splot 'test.dat' with dots; pause 5" | gnuplot

do_exponential_2d:
	./test_exponential 10000 2 > test.dat
	echo "plot 'test.dat' with dots; pause 3" | gnuplot

do_exponential_3d:
	./test_exponential 10002 3 > test.dat
	echo "set parametric; splot 'test.dat' with dots; pause 3" | gnuplot

do_sobol_2d:
	./test_sobol 10000 2 > test.dat
	echo "plot 'test.dat' with dots; pause 3" | gnuplot

do_sobol_3d:
	./test_sobol 10002 3 > test.dat
	echo "set parametric; splot 'test.dat' with dots; pause 3" | gnuplot

do_halton_2d:
	./test_halton 10000 2 > test.dat
	echo "plot 'test.dat' with dots; pause 3" | gnuplot

do_halton_3d:
	./test_halton 10002 3 > test.dat
	echo "set parametric; splot 'test.dat' with dots; pause 3" | gnuplot

do_hammersley_2d:
	./test_hammersley 10000 2 > test.dat
	echo "plot 'test.dat' with dots; pause 3" | gnuplot

do_hammersley_3d:
	./test_hammersley 10002 3 > test.dat
	echo "set parametric; splot 'test.dat' with dots; pause 3" | gnuplot

do_jittered_halton_2d:
	./test_jittered_halton 10000 2 > test.dat
	echo "plot 'test.dat' with dots; pause 3" | gnuplot

do_jittered_halton_3d:
	./test_jittered_halton 10002 3 > test.dat
	echo "set parametric; splot 'test.dat' with dots; pause 3" | gnuplot

do_jittered_hammersley_2d:
	./test_jittered_hammersley 10000 2 > test.dat
	echo "plot 'test.dat' with dots; pause 3" | gnuplot

do_jittered_hammersley_3d:
	./test_jittered_hammersley 10002 3 > test.dat
	echo "set parametric; splot 'test.dat' with dots; pause 3" | gnuplot

clean: 
	-rm -f *.o $(PROGS)

depend:
	$(CXX) -MM $(SOURCES) > .depend
