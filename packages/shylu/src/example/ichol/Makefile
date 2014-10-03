info :
	@echo "make one TEST=test"
	@echo " - where test is one of $(TESTS)"
	@echo "make all"
	@echo "make run"

OBJ   = 

TEST  = ichol_left_blocked
TESTS = \
	crs_matrix_base \
	crs_matrix_view \
	crs_row_view \
	partition \
	dot \
	ichol_left_unblocked \
	ichol_left_blocked \
	scotch

TRILINOS_DIR = ${HOME}/Work/lib/trilinos/kokkos-dev
QTHREADS_DIR = ${HOME}/Work/lib/qthreads
SCOTCH_DIR   = ${HOME}/Work/lib/scotch/6.0.0

include $(TRILINOS_DIR)/include/Makefile.export.Kokkos
include $(TRILINOS_DIR)/include/Makefile.export.KokkosExample

CC = g++
CCFLAGS = -g -O0 -std=c++11

CCINCS = -I. -I${QTHREADS_DIR}/include -I${TRILINOS_DIR}/include -I${SCOTCH_DIR}/include
CCWARN = -Wno-pointer-arith
CCDEFS = -D KOKKOS_DONT_INCLUDE_CORE_CONFIG_H \
	 -D KOKKOS_HAVE_PTHREAD \
	 -D KOKKOS_HAVE_OPENMP 

CCDEFS =

LDFLAGS = -L$(QTHREADS_DIR)/lib \
	  -lqthread \
	  -L${SCOTCH_DIR}/lib \
	  -lscotch -lscotcherr -lscotcherrexit \
	  -L${TRILINOS_DIR}/lib \
          -lkokkoscorelinalg \
	  -lkokkoscontainers \
	  -lkokkoscore \
	  -lkokkos \
	  -lkokkosTPL_unused_dummy \
	  ${KokkosExample_TPL_LIBRARIES} \
	  -lpthread

CPP_WORK = g++ $(CCFLAGS) $(CCINCS) $(CCDEBUG) $(CCWARN) $(CCDEFS)

%.o : %.cu
	@echo "Compiling $<"
	$(CPP_WORK) -o $@ -c $<

%.o : %.cpp
	@echo "Compiling $<"
	$(CPP_WORK) -o $@ -c $<

%.d : %.cpp
	@echo "Depending $<"
	$(CPP_WORK) -M $< > $@

one : $(TEST).o ${OBJ}
	@echo "Linking   $<"
	$(CPP_WORK) -o $(TEST) $(TEST).o ${OBJ} $(LDFLAGS)

all :
	@for e in $(TESTS) ; do make one TEST=$$e ; done

clean :
	@/bin/rm -f *.o *~ *.log
	@for t in $(TESTS) ; do rm -f $$t ; done
	@for d in $(DIRS); do (cd $$d; make clean ); done


DEPENDS = $(OBJ:.o=.d)
sinclude $(DEPENDS)
