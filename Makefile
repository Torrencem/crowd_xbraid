
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
   MPICXX = mpiCC
   MPIF90 = mpif90
   LFLAGS = -lm
   EXTRAFLAGS = /usr/lib/x86_64-linux-gnu/liblapacke.a /usr/lib/x86_64-linux-gnu/liblapack.a
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
   MPICXX = mpiCC
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

.PHONY: all xbraid model_problem crowd clean

.SUFFIXES:
.SUFFIXES: .c .o

all: xbraid crowd model_problem model_problem_upwind burger

xbraid: ./xbraid/braid/*.c
	cd xbraid; $(MAKE) braid

crowd: src/crowd.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) -L. -llapacke -llapack -lopenblas -lgfortran $(BRAID_FLAGS) -o crowd src/crowd.c $(BRAID_LIB_FILE) $(LFLAGS) $(EXTRAFLAGS)

utils: src/utils.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) -L. -llapacke -llapack -lopenblas -lgfortran $(BRAID_FLAGS) -o utils src/utils.c $(BRAID_LIB_FILE) $(LFLAGS) $(EXTRAFLAGS)

test: src/utils.c $(BRAID_LIB_FILE)
	@echo "Running tests in utils.c..."
	$(MPICC) -D TESTS -D NOINLINE $(CFLAGS) -Wextra -L. -llapacke -llapack -lopenblas -lgfortran $(BRAID_FLAGS) -o tests src/utils.c $(BRAID_LIB_FILE) $(LFLAGS) $(EXTRAFLAGS) -g -O0
	./tests

testxx: src/utils.cpp $(BRAID_LIB_FILE)
	@echo "Running tests in utils.cpp..."
	$(MPICXX) -D TESTS -D NOINLINE $(CXXFLAGS) -Wextra -L. -llapacke -llapack -lopenblas -lgfortran $(BRAID_FLAGS) -o tests src/utils.cpp $(BRAID_LIB_FILE) $(LFLAGS) $(EXTRAFLAGS) -g -O0
	./tests

burger: src/burger.cpp $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICXX) $(CXXFLAGS) -Wextra -L. -llapacke -llapack -lopenblas -lgfortran $(BRAID_FLAGS) -o burger src/burger.cpp $(BRAID_LIB_FILE) $(LFLAGS) $(EXTRAFLAGS)

model_problem: src/model_problem.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) -L. -llapacke -llapack -lopenblas -lgfortran $(BRAID_FLAGS) -o model_problem src/model_problem.c $(BRAID_LIB_FILE) $(LFLAGS) $(EXTRAFLAGS)

model_problem_backwards: src/model_problem.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) -D MODEL_BACKWARDS -L. -llapacke -llapack -lopenblas -lgfortran $(BRAID_FLAGS) -o model_problem_backwards src/model_problem.c $(BRAID_LIB_FILE) $(LFLAGS) $(EXTRAFLAGS)

model_problem_upwind: src/model_problem_upwind.c $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICC) $(CFLAGS) -L. -llapacke -llapack -lopenblas -lgfortran $(BRAID_FLAGS) -o model_problem_upwind src/model_problem_upwind.c $(BRAID_LIB_FILE) $(LFLAGS) $(EXTRAFLAGS)

model_problem_pp: src/model_problem.cpp $(BRAID_LIB_FILE)
	@echo "Building" $@ "..."
	$(MPICXX) $(CXXFLAGS) -L. -llapacke -llapack -lopenblas -lgfortran $(BRAID_FLAGS) -o model_problem_pp src/model_problem.cpp $(BRAID_LIB_FILE) $(LFLAGS) $(EXTRAFLAGS)

crowd_horesh_xbraid: src/crowd_horesh_xbraid.cpp
	@echo "Building" $@ "..."
	$(MPICXX) $(CXXFLAGS) -L. -llapacke -llapack -lopenblas -lgfortran $(BRAID_FLAGS) -o crowd_horesh_xbraid src/crowd_horesh_xbraid.cpp -I $(EIGEN_DIR) $(BRAID_LIB_FILE) $(LFLAGS) $(EXTRAFLAGS)

crowd_horesh: src/crowd_horesh.cpp
	@echo "Building" $@ "..."
	clang++ $(CXXFLAGS) -o crowd_horesh src/crowd_horesh.cpp -I $(EIGEN_DIR) $(EXTRAFLAGS) -Wextra -std=c++14

clean:
	rm -f *.out.*
	rm -f *.o crowd utils model_problem tests model_problem_backwards model_problem_upwind burger

fmt: src/*.c
	for file in $^ ; do \
		echo "Formatting" $${file} "..." ; \
		clang-format $${file} > $${file}.fmt ; \
		mv $${file}.fmt $${file} ; \
	done

cleanall:
	$(MAKE) clean
	cd xbraid; $(MAKE) clean
