#
# a basic makefile to build examples
#

# Location of zoltan library
BLD=../../../../build/packages/zoltan/src

# Location of zoltan header files
SRC=../../src

# Preferred compiler
CC=mpicc

# MPI include and library options if needed
MPI_INCLUDE=
MPI_LIBS=


CFLAGS=-g -Wall

all:   simpleBLOCK  simpleRCB


%: %.c
	$(CC) $(CFLAGS) $(MPI_INCLUDE) -I $(SRC)/include -o $@ $< $(BLD)/libzoltan.a -lm $(MPI_LIBS)

clean:
	@rm -rf simpleBLOCK simpleRCB
