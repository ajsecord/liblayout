COMMON_CFLAGS=-g -Wall -Wno-unused-function
COMMON_CPPFLAGS=-I. -Iext 

# Debug, non-optimized
CFLAGS=$(COMMON_CFLAGS)
CPPFLAGS=$(COMMON_CPPFLAGS)

# Release, optimized
#CFLAGS=$(COMMON_CFLAGS) -Os
#CPPFLAGS=$(COMMON_CPPFLAGS) -DNDEBUG

VPATH=src ext/macopt/newansi ext/random

all: test

test: test.o random.o liblayout.a libmacopt.a
	$(CC) -o $@ $? -framework OpenGL -framework GLUT

liblayout.a: liblayout.a(layout.o overlap.o)
	ranlib $@

libmacopt.a: libmacopt.a(macopt.o nrutil.o r.o)
	ranlib $@

clean: 
	-rm -f test liblayout.a libmacopt.a *.o

docs:
	doxygen doc/Doxygen

realclean: clean
	-rm -rf doc/code/html doc/code/latex
