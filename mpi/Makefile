# your choice of compiler
CC = mpicc

# Add your choice of flags
CFLAGS = -O3 -Wall -Wextra -fopt-info-vec -g -fopenmp-simd
LDLIBS = -lm -fopenmp-simd

all : cg

cg : cg.o mmio.o
mmio.o : mmio.c mmio.h
cg.o : cg.c mmio.h

.PHONY: clean
clean :
	rm -rf *.o cg

git :
	git add -A *
	git commit -a -m"$m"
	git push -u origin master
