#
# a simple makefile to build comm_main.c
#
# To run it:
#
# ./comm_main test/comm_input.dat
#

BLD=../../../../../build/packages/zoltan/src
SRC=../../../src
UTILSRC=$(SRC)/Utilities

all:
	mpicc \
      -I $(SRC)/include \
      -I $(UTILSRC)/shared \
      -I $(UTILSRC)/Communication \
      -I $(BLD) \
      -o comm_main comm_main.c \
      $(BLD)/libzoltan.a
