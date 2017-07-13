#Settings for what to build
KOKKOS_DEVICES = "Pthread"
KOKKOS_ARCH = "SNB,Kepler35"

# To enable Trilinos you need to copy KokkosKernels_config.h and KokkosKernels_ETIHelperMacros.h 
# from an appropriately configure Trilinos build directory (you don't need to build).
# The rest is accessed from Trilinos source.
BUILD_TRILINOS = 0
BUILD_CUSPARSE = 0
BUILD_MKL = 0

#Path settings
TRILINOS_PATH = ${HOME}/Trilinos
KOKKOS_PATH = ${HOME}/Kokkos/kokkos
KOKKOS_KERNELS_PATH = ${HOME}/Kokkos/kokkos-kernels

#MKL_PATH=/home/projects/install/rhel6-x86_64/sems/compiler/intel/16.0.1/base/mkl
#On Sandia Testbeds this should work:
#MKL_PATH=${MKLROOT}
#On Sandia SEMS modules this should work:
MKL_PATH=${SEMS_INTEL_ROOT}/mkl
MKL_LIBS=-L${MKL_PATH}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core

SRC = $(wildcard ${KOKKOS_KERNELS_PATH}/example/fenl/*.cpp)
HEADERS = $(wildcard ${KOKKOS_KERNELS_PATH}/example/fenl/*.hpp) 

tmp := $(shell echo "\#include<KokkosKernels_ETIMacros.hpp>" > KokkosKernels_ETIHelperMacros.tmp )

KOKKOS_INTERNAL_LS_CONFIG := $(shell ls KokkosKernels_ETIHelperMacros.h)
ifeq ($(KOKKOS_INTERNAL_LS_CONFIG), KokkosKernels_ETIHelperMacros.h)
	KOKKOS_INTERNAL_NEW_CONFIG := $(strip $(shell diff KokkosKernels_ETIHelperMacros.h KokkosKernels_ETIHelperMacros.tmp | grep define | wc -l))
else
	KOKKOS_INTERNAL_NEW_CONFIG := 1
endif

ifneq ($(KOKKOS_INTERNAL_NEW_CONFIG), 0)
	tmp := $(shell cp KokkosKernels_ETIHelperMacros.tmp KokkosKernels_ETIHelperMacros.h)
	tmp := $(shell touch KokkosKernels_config.h)
endif

default: build
	echo "Start Build"

CXXFLAGS = -DMAKE_BUILD -g -O3 -I./converter -I./ -I${KOKKOS_KERNELS_PATH}/example/fenl -I${KOKKOS_KERNELS_PATH}/src -I${KOKKOS_KERNELS_PATH}/src/impl
LINKFLAGS = -O3
LIB = 

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/config/nvcc_wrapper
EXE = fenl.cuda
else
CXX = g++
EXE = fenl.host
endif

LINK = ${CXX}

ifeq (${BUILD_MKL}, 1)
CXXFLAGS += -I${MKL_PATH}/include -DHAVE_MKL -mkl
LIB += ${MKL_LIBS} -mkl
endif

ifeq (${BUILD_TRILINOS}, 1)
CXXFLAGS += -DHAVE_TRILINOS
SRC += $(wildcard ${KOKKOS_KERNELS_PATH}/src/impl/*spmv*.cpp)
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

%.o:%.cpp $(KOKKOS_CPP_DEPENDS) KokkosKernels_config.h KokkosKernels_ETIHelperMacros.h ${HEADERS} 
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $<
