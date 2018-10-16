KOKKOS_PATH = ${HOME}/Kokkos/kokkos
KOKKOS_SRC_PATH = ${KOKKOS_PATH}
KOKKOSKERNELS_PATH = ${HOME}/Kokkos/kokkos-kernels
KOKKOSKERNELS_SRC_PATH = ${KOKKOSKERNELS_PATH}
GTEST_PATH = ${KOKKOS_SRC_PATH}/tpls/gtest

KOKKOS_DEVICES = "OpenMP"
KOKKOS_ARCH = "SNB,Volta70"

KOKKOSKERNELS_SCALARS =

SRC = $(wildcard *.cpp)
HDR = $(wildcard *.hpp)
KOKKOS_DEVICES=OpenMP
KOKKOS_CUDA_OPTIONS=enable_lambda

default: build
	echo "Start Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
else
CXX = g++
endif

CXXFLAGS = -O3 -g
LINK = ${CXX}
LINKFLAGS = 
EXE = test.exe 

DEPFLAGS = -M

OBJ = $(SRC:.cpp=.o)
LIB =

OBJ += gtest-all.o

INC  = -I${GTEST_PATH} -I$(KOKKOSKERNELS_SRC_PATH)/unit_test/
INC += -I$(KOKKOSKERNELS_SRC_PATH)/unit_test/blas
INC += -I$(KOKKOSKERNELS_SRC_PATH)/unit_test/sparse
INC += -I$(KOKKOSKERNELS_SRC_PATH)/unit_test/graph
INC += -I$(KOKKOSKERNELS_SRC_PATH)/unit_test/../test_common
INC += -I$(KOKKOSKERNELS_SRC_PATH)/unit_test/batched
INC += -I$(KOKKOSKERNELS_SRC_PATH)/unit_test/cuda
INC += -I$(KOKKOSKERNELS_SRC_PATH)/unit_test/serial
INC += -I$(KOKKOSKERNELS_SRC_PATH)/unit_test/threads
INC += -I$(KOKKOSKERNELS_SRC_PATH)/unit_test/openmp

TEST_HEADERS =  $(wildcard $(KOKKOSKERNELS_SRC_PATH)/unit_test/*.hpp)
TEST_HEADERS += $(wildcard $(KOKKOSKERNELS_SRC_PATH)/unit_test/blas/*.hpp)
TEST_HEADERS += $(wildcard $(KOKKOSKERNELS_SRC_PATH)/unit_test/sparse/*.hpp)
TEST_HEADERS += $(wildcard $(KOKKOSKERNELS_SRC_PATH)/unit_test/graph/*.hpp)
TEST_HEADERS += $(wildcard $(KOKKOSKERNELS_SRC_PATH)/unit_test/../test_common/*.hpp)
TEST_HEADERS += $(wildcard $(KOKKOSKERNELS_SRC_PATH)/unit_test/batched/*.hpp)

include $(KOKKOS_PATH)/Makefile.kokkos
include ${KOKKOSKERNELS_PATH}/Makefile.kokkos-kernels

build: $(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS) $(KOKKOSKERNELS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(KOKKOSKERNELS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(KOKKOSKERNELS_LIBS) $(LIB) -o $(EXE)


clean: kokkos-clean 
	rm -f *.o *.cuda *.host

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS) $(HDR)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(KOKKOSKERNELS_CPPFLAGS) $(INC) -c $<

gtest-all.o:$(GTEST_PATH)/gtest/gtest-all.cc
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(INC) $(CXXFLAGS) -c $(GTEST_PATH)/gtest/gtest-all.cc
