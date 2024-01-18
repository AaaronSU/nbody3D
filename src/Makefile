# clang -march=native -g3 -Ofast -fopenmp=libgomp nbody.c -o nbody3D -lm -I/usr/lib/gcc/x86_64-linux-gnu/11/include -L/OPENMP_RUNTIME_LIB

CC=gcc

CFLAGS=-march=native -g3
# -fopt-info-all=nbody.gcc.optrpt
OFLAGS=-Ofast -fopenmp 

all: nbody3D
	taskset -c $(NP) ./nbody
	taskset -c $(NP) ./nbody2
	taskset -c $(NP) ./nbody3

nbody3D: nbody.c nbody2.c
	$(CC) $(CFLAGS) $(OFLAGS) nbody.c -o nbody -lm
	$(CC) $(CFLAGS) $(OFLAGS) nbody2.c -o nbody2 -lm
	$(CC) $(CFLAGS) $(OFLAGS) nbody3.c -o nbody3 -lm

clean:
	rm -Rf *~ nbody3D *.optrpt nbody nbody2 nbody3