#!/bin/bash

# CLANG_CFLAGS
CLANG_CFLAGS="-march=native -g3 -fopenmp=libgomp -lm -I/usr/lib/gcc/x86_64-linux-gnu/11/include -L/OPENMP_RUNTIME_LIB"

# GCC_CFLAGS
GCC_CFLAGS="-march=native -g3 -fopenmp -funroll-loops"

echo "gcc_o2"
gcc $GCC_CFLAGS -O2 nbody5.c -o bin/gcc_O2_nbody -lm && bin/gcc_O2_nbody > data/gcc_O2_nbody
cat data/gcc_O2_nbody

echo "gcc_o3"
gcc $GCC_CFLAGS -O3 nbody5.c -o bin/gcc_O3_nbody -lm && bin/gcc_O3_nbody > data/gcc_O3_nbody 
cat data/gcc_O3_nbody

echo "gcc_ofast"
gcc $GCC_CFLAGS -Ofast nbody5.c -o bin/gcc_Ofast_nbody -lm && bin/gcc_Ofast_nbody > data/gcc_Ofast_nbody
cat data/gcc_Ofast_nbody

echo "clang_o2"
clang $CLANG_CFLAGS -O2 nbody5.c -o bin/clang_O2_nbody -lm && bin/clang_O2_nbody > data/clang_O2_nbody
cat data/clang_O2_nbody

echo "clang_o3"
clang $CLANG_CFLAGS -O3 nbody5.c -o bin/clang_O3_nbody -lm && bin/clang_O3_nbody > data/clang_O3_nbody
cat data/clang_O3_nbody

echo "clang_ofast"
clang $CLANG_CFLAGS -Ofast nbody5.c -o bin/clang_Ofast_nbody -lm && bin/clang_Ofast_nbody > data/clang_Ofast_nbody
cat data/clang_Ofast_nbody