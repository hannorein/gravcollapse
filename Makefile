PROBLEMDIR=$(shell basename `dirname \`pwd\``)"/"$(shell basename `pwd`)
export OPT=-O3 
export OPENGL=1
export LIBPNG=0
export MPI=0
export OPENMP=1
export CC=gcc
export OMP_NUM_THREADS = 8 * 20
CFLAGS=-Wall -g

all:
	# Setup link to different modules
	ln -fs gravity_tree.c ../../src/gravity.c
	ln -fs boundaries_shear.c ../../src/boundaries.c
	ln -fs collisions_tree.c ../../src/collisions.c
	# Setup link to problem file
	ln -fs ../$(PROBLEMDIR)/gravcollapse.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	# Copy result
	cp ../../src/rebound .

clean:
	$(MAKE) -C ../../src/ clean
	rm -vf rebound
