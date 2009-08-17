#
# a basic makefile to build examples
#

# Location of zoltan library, also location of Zoltan_config.h if
# using the autoconf build harness

BLD=../../../../build/packages/zoltan/src

# Location of zoltan "src" directory

SRC=../../src

# Location of zoltan header files

ZOLTAN_INCLUDE=-I$(SRC)/include

# Headers for zoltan utilities, not often used by applications
# Needed for Zoltan_Hash in simpleGRAPH.c

ZOLTAN_UTIL_INCLUDE=-I$(SRC)/zz  \
                    -I$(SRC)/par \
                    -I$(SRC)/order \
                    -I$(SRC)/lb \
                    -I$(SRC)/params \
                    -I$(SRC)/Utilities/shared

# Preferred compiler
CC=mpicc

# MPI include and library options if needed
MPI_INCLUDE=
MPI_LIBS=


CFLAGS=-g -Wall

all:   simpleBLOCK  simpleRCB simpleGRAPH


%: %.c
	$(CC) $(CFLAGS) -I$(BLD) $(MPI_INCLUDE) $(ZOLTAN_INCLUDE) $(ZOLTAN_UTIL_INCLUDE) \
        -o $@ $< $(BLD)/libzoltan.a -lm $(MPI_LIBS)

clean:
	@rm -rf $(all)
