
#=======================================================================
#=================== Settings ==========================================
#=======================================================================

KOKKOS_PATH ?= ${HOME}/Kokkos/kokkos
KOKKOSKERNELS_PATH ?= ${HOME}/Kokkos/kokkos-kernels
KOKKOSKERNELS_SRC_PATH ?= ${KOKKOSKERNELS_PATH}
GTEST_PATH = ${KOKKOS_PATH}/tpls/gtest

KOKKOS_DEVICES=OpenMP
KOKKOS_CUDA_OPTIONS=enable_lambda

KOKKOSKERNELS_SCALARS="'double,complex<double>'"
KOKKOSKERNELS_OPTIONS=eti-only


CXX = g++

CXXFLAGS = -O3 -g
LINK = ${CXX}
LINKFLAGS =  

DEPFLAGS = -M

OBJ =
LIB = -lpthread


#=======================================================================
#=================== Include Kokkos and KokkosKernels ==================
#=======================================================================

#include $(KOKKOS_PATH)/Makefile.kokkos
include $(KOKKOSKERNELS_PATH)/Makefile.kokkos-kernels

#=======================================================================
#=================== Add common source files ===========================
#=======================================================================

EXECUTABLES =
TEST_HEADERS =
TEST_TARGETS =
TARGETS =

#=======================================================================
#= Example Directories
#=======================================================================
TEST_HEADERS += $(wildcard $(KOKKOSKERNELS_SRC_PATH)/example/hashmap_accumulator/*.hpp)
EXECUTABLES += $(wildcard ${KOKKOSKERNELS_SRC_PATH}/example/hashmap_accumulator/*.cpp) 

#EXECUTABLES += $(wildcard ${KOKKOSKERNELS_SRC_PATH}/example/graph/PartitioningExample.cpp) 

#=======================================================================
#========================== TARGETS ====================================
#=======================================================================

vpath %.cpp $(sort $(dir $(EXECUTABLES)))
OBJS = $(notdir $(EXECUTABLES:.cpp=.o))
TEST_TARGETS = $(notdir $(EXECUTABLES:.cpp=.exe))


#=======================================================================
#=================== Make Rules ========================================
#=======================================================================


clean: kokkos-clean kokkoskernels-clean
	rm -f *.o *.exe

default: $(TEST_TARGETS)

build: $(TEST_TARGETS)

%.exe:%.o $(KOKKOS_LINK_DEPENDS) $(KOKKOSKERNELS_LINK_DEPENDS) $(TEST_HEADERS)
	$(LINK) $(EXTRA_PATH) $< $(KOKKOSKERNELS_LDFLAGS) $(KOKKOSKERNELS_LIBS) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o $@

%.o:%.cpp $(KOKKOS_CPP_DEPENDS) $(KOKKOSKERNELS_CPP_DEPENDS) $(TEST_HEADERS) 
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOSKERNELS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(INC) $(CXXFLAGS) $(EXTRA_INC) -I. -c $< -o $(notdir $@)

#depend:
#	makedepend -Y ${EXECUTABLES} $(TEST_HEADERS)
# DO NOT DELETE
