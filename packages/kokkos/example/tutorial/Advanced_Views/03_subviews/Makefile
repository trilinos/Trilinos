KOKKOS_PATH = ../../../..
SRC = $(wildcard *.cpp)

ifeq ($(CUDA), yes)
CC = nvcc
CCFLAGS = -O3 -arch=sm_35  
LINK = g++
LINKFLAGS = 
EXE = $(SRC:.cpp=.cuda)
else
CC = g++
CCFLAGS = -O3
LINK = g++
LINKFLAGS =  
EXE = $(SRC:.cpp=.host)
endif

DEPFLAGS = -M

OBJ = $(SRC:.cpp=.o)
LIB =

include $(KOKKOS_PATH)/Makefile.kokkos

EXTRA_INC = $(KOKKOS_INC)
EXTRA_LIB = $(KOKKOS_LINK)

all: $(OBJ)
	$(LINK) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(EXTRA_LIB) $(LIB) -o $(EXE)
	rm *.o *.d libkokkoscore.a kokkos_depend.cpp

clean: 
	rm *.o *.d libkokkoscore.a kokkos_depend.cpp *.cuda *.host

# Compilation rules

%.o:%.cu
	$(CC) $(CCFLAGS) $(EXTRA_INC) -c $<

%.o:%.cpp
	$(CC) $(CCFLAGS) $(EXTRA_INC) -c $<

%.d:%.cpp
	$(CC) $(CCFLAGS) $(EXTRA_INC) $(DEPFLAGS) $< > $@

# Individual dependencies

DEPENDS = $(OBJ:.o=.d)
sinclude $(DEPENDS)

