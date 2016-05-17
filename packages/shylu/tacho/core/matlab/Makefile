info :
	@echo "make one TEST=test"
	@echo " - where test is one of $(TESTS)"
	@echo "make all"
	@echo "make run"

SRC   = 

OBJ   = $(SRC:.cpp=.o)

TEST  = ichol_strong_scale
TESTS = \
	crs_matrix_base \
	crs_matrix_view \
	crs_row_view \
	\
	scotch \
	qthread \
	crs_matrix_helper \
	\
	partition \
	dot \
	\
	ichol_left_unblocked \
	ichol_left_blocked \
	\
	ichol_right_unblocked \
	ichol_right_blocked \
	\
	crs_matrix_base_future \
	task_policy_graphviz \
	crs_task_view \
	\
	ichol_left_by_blocks \
	ichol_right_by_blocks \
	\
	ichol_strong_scale

MEXFILE  = mex/ichol_unblocked
MEXFILES = \
	mex/kokkos \
	mex/ichol_unblocked \
	mex/ichol_blocked 

TRILINOS_DIR = ${HOME}/Work/lib/trilinos/kokkos-openmp
QTHREADS_DIR = ${HOME}/Work/lib/qthreads
SCOTCH_DIR   = ${HOME}/Work/lib/scotch/6.0.0

include $(TRILINOS_DIR)/include/Makefile.export.Kokkos
include $(TRILINOS_DIR)/include/Makefile.export.KokkosExample

CC = g++
MEX = mex

CCFLAGS = -g -O3 -std=c++11 -fopenmp -fPIC -fdiagnostics-color=auto 
#-march=native 

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

CPP_WORK = $(CC) $(CCFLAGS) $(CCINCS) $(CCDEBUG) $(CCWARN) $(CCDEFS)
MEX_WORK = $(MEX) -v -largeArrayDims \
	CC="${CC}" \
	CXXFLAGS="${CCFLAGS}" \
	CXXOPTIMFLAGS="${CCFLAGS} -DNDEBUG" \
	LD="\${CC} -fopenmp -O3" \
	$(CCINCS) $(CCDEBUG) $(CCDEFS)

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

matlab : $(MEXFILE).cpp
	@echo "Linking   $<"
	$(MEX_WORK) $(MEXFILE).cpp ${OBJ} $(LDFLAGS)

matlab-all : 
	@for e in $(MEXFILES) ; do make matlab MEXFILE=$$e ; done

clean :
	@/bin/rm -f *.o *~ *.log *.mexa64 
	@for t in $(TESTS) ; do rm -f $$t ; done


DEPENDS = $(OBJ:.o=.d) 
sinclude $(DEPENDS)
