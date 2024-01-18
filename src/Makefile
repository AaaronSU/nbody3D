# clang -march=native -g3 -Ofast -fopenmp=libgomp nbody.c -o nbody3D -lm -I/usr/lib/gcc/x86_64-linux-gnu/11/include -L/OPENMP_RUNTIME_LIB

CC=gcc
NP=3

CFLAGS=-march=native -g3
# -fopt-info-all=nbody.gcc.optrpt
OFLAGS=-O3 -fopenmp -fopt-info-all=nbody.gcc.optrpt

all: nbody3D
	# taskset -c $(NP) ./bin/nbody
	# taskset -c $(NP) ./bin/nbody2
	# taskset -c $(NP) ./bin/nbody3
	# taskset -c $(NP) ./bin/nbody4
	taskset -c $(NP) ./bin/nbody5

nbody3D: nbody.c nbody2.c nbody3.c nbody4.c nbody5.c
	$(CC) $(CFLAGS) $(OFLAGS) nbody.c -o bin/nbody -lm
	$(CC) $(CFLAGS) $(OFLAGS) nbody2.c -o bin/nbody2 -lm
	$(CC) $(CFLAGS) $(OFLAGS) nbody3.c -o bin/nbody3 -lm
	$(CC) $(CFLAGS) $(OFLAGS) nbody4.c -o bin/nbody4 -lm
	$(CC) $(CFLAGS) $(OFLAGS) nbody5.c -o bin/nbody5 -lm

clean:
	rm -Rf *~ nbody3D *.optrpt bin/*