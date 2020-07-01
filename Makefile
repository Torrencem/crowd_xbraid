
# Imported from makefile.inc in XBraid repository
HOSTNAME := $(shell hostname)

ifeq ($(shell uname -s), Darwin)
   # for Jacob's MacBook running Homebrew 
   # Need to specifically include lstdc++ (!!)
   CC = gcc
   MPICC = mpicc
   MPICXX = mpicxx
   MPIF90 = mpif90
   LFLAGS = -lm -lstdc++
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall -L/usr/local/opt/openblas/lib -I/usr/local/opt/openblas/include -L/usr/local/opt/lapack/lib -I/usr/local/opt/lapack/insude
      CXXFLAGS = -g -Wall
      FORTFLAGS = -Og -Wall
   else
      CFLAGS = -O -Wall -L/usr/local/opt/openblas/lib -I/usr/local/opt/openblas/include -L/usr/local/opt/lapack/lib -I/usr/local/opt/lapack/insude
      CXXFLAGS = -O -Wall
      FORTFLAGS = -O1 -Wall
   endif
else ifeq ($(findstring cab,$(HOSTNAME)),cab)
   # for Cab
   MPICC = mpiicc
   MPICXX = mpiicpc
   MPIF90 = mpif90
   LFLAGS = -lm
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall
      CXXFLAGS = -g -Wall
      FORTFLAGS = -Og -Wall
   else
      CFLAGS = -O -Wall
      CXXFLAGS = -O -Wall
      FORTFLAGS = -O1 -Wall
   endif   
else ifeq ($(findstring vulcan,$(HOSTNAME)),vulcan)
   # for Vulcan
   MPICC = mpixlc
   MPICXX = mpixlcxx
   MPIF90 = mpixlf90
   LFLAGS = -lm
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall
      CXXFLAGS = -g -Wall
      FORTFLAGS = -Og -Wall
   else
      CFLAGS = -O -Wall
      CXXFLAGS = -O -Wall
      FORTFLAGS = -O1 -Wall
   endif
else ifeq ($(shell uname -s),Linux)
   MPICC = mpicc
   MPICXX = mpiCC
   MPIF90 = mpif90
   LFLAGS = -lm
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall
      CXXFLAGS = -g -Wall
      FORTFLAGS = -g -Wall
   else
      CFLAGS = -O -Wall
      CXXFLAGS = -O -Wall
      FORTFLAGS = -O1 -Wall
   endif
else
   MPICC = mpicc
   MPICXX = mpiCC
   MPIF90 = mpif90
   LFLAGS = -lm
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall
      CXXFLAGS = -g -Wall
      FORTFLAGS = -g -Wall
   else
      CFLAGS = -O -Wall
      CXXFLAGS = -O -Wall
      FORTFLAGS = -O1 -Wall
   endif
endif

# Modified from examples Makefile in XBraid repository
BRAID_DIR = ./xbraid/braid
BRAID_FLAGS = -I$(BRAID_DIR)
BRAID_LIB_FILE = $(BRAID_DIR)/libbraid.a

.PHONY: all xbraid crowd clean

.SUFFIXES:
.SUFFIXES: .c .o

all: xbraid crowd

xbraid: ./xbraid/braid/*.c
	cd xbraid; $(MAKE) braid

crowd: src/main.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) -L. -llapacke -llapack -lopenblas -lgfortran $(BRAID_FLAGS) -o crowd src/main.c $(BRAID_LIB_FILE) $(LFLAGS)

clean:
	rm -f *.out.*
	rm -f *.o crowd
	cd xbraid; $(MAKE) clean

cleanall:
	$(MAKE) clean
	cd lapack; $(MAKE) clean
