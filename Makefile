
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
   EXTRAFLAGS = 
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall -L/usr/local/opt/openblas/lib -I/usr/local/opt/openblas/include -L/usr/local/opt/lapack/lib -I/usr/local/opt/lapack/insude -D NOINLINE
      CXXFLAGS = -g -Wall -L/usr/local/opt/openblas/lib -I/usr/local/opt/openblas/include -L/usr/local/opt/lapack/lib -I/usr/local/opt/lapack/insude -D NOINLINE
      FORTFLAGS = -Og -Wall
   else
      CFLAGS = -O -Wall -L/usr/local/opt/openblas/lib -I/usr/local/opt/openblas/include -L/usr/local/opt/lapack/lib -I/usr/local/opt/lapack/insude
      CXXFLAGS = -O -Wall -L/usr/local/opt/openblas/lib -I/usr/local/opt/openblas/include -L/usr/local/opt/lapack/lib -I/usr/local/opt/lapack/insude
      FORTFLAGS = -O1 -Wall
   endif
else ifeq ($(findstring cab,$(HOSTNAME)),cab)
   # for Cab
   MPICC = mpiicc
   MPICXX = mpiicpc
   MPIF90 = mpif90
   LFLAGS = -lm
   EXTRAFLAGS = 
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall -D NOINLINE
      CXXFLAGS = -g -Wall -D NOINLINE
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
   EXTRAFLAGS = 
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall -D NOINLINE
      CXXFLAGS = -g -Wall -D NOINLINE
      FORTFLAGS = -Og -Wall
   else
      CFLAGS = -O -Wall
      CXXFLAGS = -O -Wall
      FORTFLAGS = -O1 -Wall
   endif
else ifeq ($(shell uname -s),Linux)
   MPICC = mpicc
   MPICXX = mpicxx
   MPIF90 = mpif90
   LFLAGS = -lm
   EXTRAFLAGS = -std=c++14
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall -D NOINLINE
      CXXFLAGS = -g -Wall -D NOINLINE
      FORTFLAGS = -g -Wall
   else
      CFLAGS = -O -Wall
      CXXFLAGS = -O -Wall
      FORTFLAGS = -O1 -Wall
   endif
else
   MPICC = mpicc
   MPICXX = mpicxx
   MPIF90 = mpif90
   LFLAGS = -lm
   EXTRAFLAGS = 
   ifeq ($(optlevel),DEBUG)
      CFLAGS = -g -Wall -D NOINLINE
      CXXFLAGS = -g -Wall -D NOINLINE
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

EIGEN_DIR = ./eigen/

.PHONY: all xbraid clean list

.SUFFIXES:
.SUFFIXES: .c .o .cpp

all: xbraid crowd_horesh crowd_horesh_xbraid

list:
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$'

xbraid: ./xbraid/braid/*.c
	cd xbraid; $(MAKE) braid

ex-01-mod: src/ex-01-mod.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) -o ex-01-mod src/ex-01-mod.c $(BRAID_LIB_FILE) $(LFLAGS) $(EXTRAFLAGS)

trimgrit-burgers: src/trimgrit-burgers.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) $(BRAID_FLAGS) -o trimgrit-burgers src/trimgrit-burgers.c $(BRAID_LIB_FILE) $(LFLAGS) $(EXTRAFLAGS)

crowd_horesh_xbraid: src/crowd_horesh_xbraid.cpp
	@echo "Building" $@ "..."
	$(MPICXX) $(CXXFLAGS) $(BRAID_FLAGS) -o crowd_horesh_xbraid src/crowd_horesh_xbraid.cpp -I $(EIGEN_DIR) $(BRAID_LIB_FILE) $(LFLAGS) $(EXTRAFLAGS)

crowd_horesh_xbraid_schur: src/crowd_horesh_xbraid_schur.cpp
	@echo "Building" $@ "..."
	$(MPICXX) $(CXXFLAGS) $(BRAID_FLAGS) -o crowd_horesh_xbraid_schur src/crowd_horesh_xbraid_schur.cpp -I $(EIGEN_DIR) $(BRAID_LIB_FILE) $(LFLAGS) $(EXTRAFLAGS)

crowd_horesh: src/crowd_horesh.cpp
	@echo "Building" $@ "..."
	$(MPICXX) $(CXXFLAGS) -o crowd_horesh src/crowd_horesh.cpp -I $(EIGEN_DIR) $(EXTRAFLAGS) -Wextra -std=c++14

run_horesh: crowd_horesh
	./crowd_horesh

run_horesh_xbraid: crowd_horesh_xbraid
	./crowd_horesh_xbraid

run_horesh_schur: crowd_horesh_xbraid_schur
	./crowd_horesh_xbraid_schur

run_ex_01_mod: ex-01-mod
	./ex-01-mod

clean:
	rm -f *.out.*
	rm -f *.o crowd_horesh crowd_horesh_xbraid crowd_horesh_xbraid_schur ex-01-mod

fmt: src/*.c src/*.cpp
	for file in $^ ; do \
		echo "Formatting" $${file} "..." ; \
		clang-format $${file} > $${file}.fmt ; \
		mv $${file}.fmt $${file} ; \
	done

cleanall:
	$(MAKE) clean
	cd xbraid; $(MAKE) clean
