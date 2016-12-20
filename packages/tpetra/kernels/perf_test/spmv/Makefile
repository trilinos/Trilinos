#Settings for what to build
KOKKOS_DEVICES = "OpenMP"
KOKKOS_ARCH = "SNB,Kepler35"

# To enable Trilinos you need to copy TpetraKernels_config.h and TpetraKernels_ETIHelperMacros.h 
# from an appropriately configure Trilinos build directory (you don't need to build).
# The rest is accessed from Trilinos source.
BUILD_TRILINOS = 0
BUILD_CUSPARSE = 0
BUILD_MKL = 0

#Path settings
TRILINOS_PATH = ${HOME}/Trilinos
KOKKOS_PATH = ${TRILINOS_PATH}/packages/kokkos

#MKL_PATH=/home/projects/install/rhel6-x86_64/sems/compiler/intel/16.0.1/base/mkl
#On Sandia Testbeds this should work:
#MKL_PATH=${MKLROOT}
#On Sandia SEMS modules this should work:
MKL_PATH=${SEMS_INTEL_ROOT}/mkl
MKL_LIBS=-L${MKL_PATH}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core

SRC = ${TRILINOS_PATH}/packages/tpetra/kernels/perf_test/spmv/test_crsmatrix.cpp 
HEADERS = $(wildcard ${TRILINOS_PATH}/packages/tpetra/kernels/perf_test/spmv/*.hpp) 

default: build
	echo "Start Build"

config_headers: 
	touch TpetraKernels_ETIHelperMacros.h
	touch TpetraKernels_config.h


CXXFLAGS = -DMAKE_BUILD -g -O3 -I./converter -I./ -I${TRILINOS_PATH}/packages/tpetra/kernels/perf_test/spmv -I${TRILINOS_PATH}/packages/tpetra/kernels/src -I${TRILINOS_PATH}/packages/tpetra/kernels/src/impl
LINKFLAGS = -O3
LIB = 

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/config/nvcc_wrapper
EXE = test_matvec.cuda
else
CXX = g++
EXE = test_matvec.host
endif

LINK = ${CXX}

ifeq (${BUILD_MKL}, 1)
CXXFLAGS += -I${MKL_PATH}/include -DHAVE_MKL -mkl
LIB += ${MKL_LIBS} -mkl
endif

ifeq (${BUILD_TRILINOS}, 1)
CXXFLAGS += -DHAVE_TRILINOS
SRC += $(wildcard ${TRILINOS_PATH}/packages/tpetra/kernels/src/impl/*spmv*.cpp)
endif


ifeq (${BUILD_CUSPARSE}, 1)
CXXFLAGS += -I${CUDA_ROOT}/include -DHAVE_CUSPARSE -DKOKKOS_USE_CUSPARSE
LIB += -lcusparse
endif

vpath %.cpp $(sort $(dir $(SRC)))

DEPFLAGS = -M

SRC_NODIR = $(notdir $(SRC))
OBJ = $(SRC_NODIR:.cpp=.o)

include $(KOKKOS_PATH)/Makefile.kokkos

build: $(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE)

clean: kokkos-clean 
	rm -f *.o *.cuda *.host TpetraKernels_config.h TpetraKernels_ETIHelperMacros.h

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS) config_headers ${HEADERS} 
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<
