#Settings for what to build
KOKKOS_DEVICES = "OpenMP"
KOKKOS_ARCH = "SNB,Kepler35"

# To enable Trilinos you need to copy KokkosKernels_config.h and KokkosKernels_ETIHelperMacros.h 
# from an appropriately configure Trilinos build directory (you don't need to build).
# The rest is accessed from Trilinos source.
BUILD_KK_KERNELS = 1
BUILD_CUSPARSE = 0
BUILD_MKL = 0

#Path settings
KOKKOS_PATH = ${HOME}/Kokkos/kokkos
KOKKOS_KERNELS_SRC_PATH = ${HOME}/Kokkos/kokkos-kernels
KOKKOS_KERNELS_INSTALL_PATH = ${HOME}/Kokkos/kokkos-kernels-install

#MKL_PATH=/home/projects/install/rhel6-x86_64/sems/compiler/intel/16.0.1/base/mkl
#On Sandia Testbeds this should work:
#MKL_PATH=${MKLROOT}
#On Sandia SEMS modules this should work:
MKL_PATH=${SEMS_INTEL_ROOT}/mkl
MKL_LIBS=-L${MKL_PATH}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core

SRC = ${KOKKOS_KERNELS_SRC_PATH}/perf_test/spmv/test_crsmatrix.cpp 
HEADERS = $(wildcard ${KOKKOS_KERNELS_SRC_PATH}/perf_test/spmv/*.hpp) 

default: build
	echo "Start Build"

config_headers: 
	touch KokkosKernels_ETIHelperMacros.h
	touch KokkosKernels_config.h


CXXFLAGS = -DMAKE_BUILD -g -O3 -I./converter -I./ -I${KOKKOS_KERNELS_SRC_PATH}/perf_test/spmv -I${KOKKOS_KERNELS_INSTALL_PATH}/include
LINKFLAGS = -O3
LIB = ${KOKKOS_KERNELS_INSTALL_PATH}/lib/libkokkos_kernels.a 

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
	rm -f *.o *.cuda *.host KokkosKernels_config.h KokkosKernels_ETIHelperMacros.h

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS) config_headers ${HEADERS} 
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<
