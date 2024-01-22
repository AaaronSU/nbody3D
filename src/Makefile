# clang -march=native -g3 -Ofast -fopenmp=libgomp nbody.c -o nbody3D -lm -I/usr/lib/gcc/x86_64-linux-gnu/11/include -L/OPENMP_RUNTIME_LIB

CC=gcc
NP=3

CFLAGS=-march=native -g3
# -fopt-info-all=nbody.gcc.optrpt
OFLAGS=-O3 -fopenmp -funroll-loops

all: nbody3D
	echo "comparaison en séquentiel"
	taskset -c $(NP) ./bin/nbody5
	echo "comparaison en parallèle"
	./bin/nbody5

nbody3D: nbody5.c
	mkdir -p bin/
	$(CC) $(CFLAGS) $(OFLAGS) nbody5.c -o bin/nbody5 -lm

clean:
	rm -Rf *~ nbody3D *.optrpt bin/* data/*
