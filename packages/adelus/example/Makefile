KOKKOS_PATH = ${HOME}/Trilinos/Trilinos.fork/packages/kokkos
KOKKOS_DEVICES = "Cuda"
EXE_NAME = "adelus_driver"

ADELUSDIR   ?= ${HOME}/Trilinos/install-trilinos/lib
TRILINOSDIR ?= ${HOME}/Trilinos/install-trilinos/lib

ADELUS_CFLAGS := -I${HOME}/Trilinos/install-trilinos/include
MPI_CFLAGS    := -I$(MPI_ROOT)/include

BLAS_LIBS    := -L$(ESSLLIBDIR64) -lesslsmp
MPI_LIBS     := -L$(MPI_ROOT)/lib -lmpi_ibm
ADELUS_LIBS  := $(ADELUSDIR)/libzadelus.a
TRILINOS_LIBS:= $(TRILINOSDIR)/libkokkoscore.a $(TRILINOSDIR)/libkokkoscontainers.a $(TRILINOSDIR)/libkokkosalgorithms.a $(TRILINOSDIR)/libkokkoskernels.a -ldl

SRC = $(wildcard *.cpp)
HDR = $(wildcard *.hpp)

default: build
	echo "Start Build"

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/bin/nvcc_wrapper
EXE = ${EXE_NAME}.adelus
TRILINOS_LIBS+= -lcudart -lcuda -lcublas
CXXFLAGS = --std=c++11 -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored -expt-extended-lambda -arch=sm_70 -rdc=true -I${CUDA_HOME}/include -fopenmp
LINKFLAGS = -arch=sm_70 -rdc=true -L${CUDA_HOME}/lib64 -L${CUDA_HOME}/lib -fopenmp
else
CXX = g++
EXE = ${EXE_NAME}.adelus
CXXFLAGS = --std=c++11 -fopenmp
LINKFLAGS = -fopenmp
endif

CXXFLAGS += -O3 -g
LINK = ${CXX}

OBJ = $(notdir $(SRC:.cpp=.o))
LIB =

OPTS=-DKKVIEW_API

$(warning $(OBJ) $(EXE) $(sort $(dir $(SRC))))

build: $(EXE)

$(EXE): $(OBJ)
	$(LINK) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(LIB) $(BLAS_LIBS) $(MPI_LIBS) $(ADELUS_LIBS) $(TRILINOS_LIBS) -o $(EXE)

clean:
	rm -f *.o *.cuda *.host *.a *.tmp *.x *.adelus

# Compilation rules

%.o:%.cpp $(HDR)
	$(CXX) $(CXXFLAGS) $(EXTRA_INC) $(MPI_CFLAGS) $(ADELUS_CFLAGS) $(OPTS) -c $< -o $(notdir $@)

test: $(EXE)
	./$(EXE)
