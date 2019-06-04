KOKKOS_PATH ?= ${HOME}/Kokkos/kokkos
KOKKOSKERNELS_PATH ?= ${HOME}/Kokkos/kokkos-kernels

# Turn of ETI
KOKKOSKERNELS_SCALARS =  

SRC ?= $(wildcard *.cpp)
HDR ?= $(wildcard *.hpp)

KOKKOS_DEVICES=OpenMP
KOKKOS_ARCH = "SNB,Kepler35"
KOKKOS_CUDA_OPTIONS=enable_lambda

default: build
	echo "Start Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/config/nvcc_wrapper
else
CXX = g++
endif

LINK = ${CXX}

CXXFLAGS = -O3 -g 
override CXXFLAGS += -I./
LINKFLAGS = 

EXE = test.x
DEPFLAGS = -M

vpath %.cpp $(sort $(dir $(SRC)))

OBJ = $(notdir $(SRC:.cpp=.o))
LIB =

include $(KOKKOS_PATH)/Makefile.kokkos
include ${KOKKOSKERNELS_PATH}/Makefile.kokkos-kernels

$(warning $(OBJ) $(EXE) $(sort $(dir $(SRC))))

build: $(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS) $(KOKKOSKERNELS_LINK_DEPENDS)
	$(LINK) $(KOKKOS_LDFLAGS) $(KOKKOSKERNELS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(KOKKOSKERNELS_LIBS) $(LIB) -o $(EXE)

clean: kokkos-clean 
	rm -f *.o *.cuda *.host

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS) $(HDR)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(KOKKOSKERNELS_CPPFLAGS) $(CXXFLAGS) $(EXTRA_INC) -c $< -o $(notdir $@)

