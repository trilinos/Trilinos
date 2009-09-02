##
## a basic makefile to build examples
##

## Location of zoltan library, also location of Zoltan_config.h if
## using the autoconf build harness

BLD=../../../../build/packages/zoltan/src

##
## Location of zoltan "src" directory (for header files)
##

SRC=../../src

# Location of zoltan header files

ZOLTAN_INCLUDE=-I$(SRC)/include

##
## Headers for zoltan utilities, not often used by applications
## Needed for use of Zoltan_Hash in simpleGRAPH.c
##

ZOLTAN_UTIL_INCLUDE=-I$(SRC)/zz  \
                    -I$(SRC)/par \
                    -I$(SRC)/order \
                    -I$(SRC)/lb \
                    -I$(SRC)/params \
                    -I$(SRC)/Utilities/shared

##
## If Zoltan was built with ParMETIS, provide the information required
## to link with the ParMETIS library.
##

PAR_LIB_DIR=-L/home/lafisk/system/parmetis/ParMetis-3.1.1
PAR_LIBS=-lparmetis -lmetis

##
## C compiler
##
CC=mpicc

##
## C++ compiler
##
CXX=mpicxx

# MPI include and library options if needed
MPI_INCLUDE=
MPI_LIBS=


CFLAGS=-g -Wall
CXXFLAGS=-g -Wall

C_EXAMPLES = simpleBLOCK  simpleRCB simpleGRAPH simplePHG 

CPP_EXAMPLES = cpp_simpleBLOCK  cpp_simpleRCB cpp_simpleGRAPH cpp_simplePHG 

all:   $(C_EXAMPLES) $(CPP_EXAMPLES)

%: %.c
	$(CC) $(CFLAGS) $(PAR_LIB_DIR)  \
        -I$(BLD) $(MPI_INCLUDE) $(ZOLTAN_INCLUDE) $(ZOLTAN_UTIL_INCLUDE) \
        -o $@ $< $(BLD)/libzoltan.a $(PAR_LIBS) -lm $(MPI_LIBS)

%: %.cpp
	$(CXX) $(CFLAGS) $(CXXFLAGS) $(PAR_LIB_DIR)  \
        -I$(BLD) $(MPI_INCLUDE) $(ZOLTAN_INCLUDE) $(ZOLTAN_UTIL_INCLUDE) \
        -o $@ $< $(BLD)/libzoltan.a $(PAR_LIBS) -lm $(MPI_LIBS)

clean:
	@rm -rf $(all)
