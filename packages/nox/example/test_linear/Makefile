# To run makefile:
#    1) set environment variable TRILINOS_ARCH to sgi, sun, tflop, or pclinux.
#       Other machines require an appropriate makefile.$(TRILINOS_ARCH) file.
#    2) Set TRILINOS_COMM to SERIAL or MPI
#    3) (Optional) Set TRILINOS_ID to make unique version for same 
#       architecture and communication mode.
#
#    4) Make the archive $(LIBAZTEC) by typing 'make'.
#

include ../../build/Makefile.defs

LIB_PATHS= $(LIBNLSPACK) $(LIBEPETRA) $(LIBAZTEC) $(LIBIFPACK) \
	$(LIBLAPACK) $(LIBSPBLAS) $(LIBBLAS) $(LIBY12M) 

#=======================================================================
# Petra test source files
#=======================================================================

TEST_CC = cc_main.cc fill.cc basis.cc NLS_Interface.cc
TEST_C = 
TEST_F = 

#=======================================================================
# TEST include files
#=======================================================================

TEST_INC = 

TEST_OBJ          = $(TEST_CC:.cc=.o) $(TEST_C:.c=.o) $(TEST_F:.f=.o)


TARGET_MPI = test_mpi
TARGET_SERIAL = test_serial
TARGET = $(TARGET_$(TRILINOS_COMM))


$(TARGET): $(TEST_OBJ) $(LIBNLSPACK)
	@echo "Making " $(TARGET)
	$(LINKER) $(ARCH_FLAGS) $(TEST_OBJ) $(LIB_PATHS) $(ARCH_LIBS) \
	$(LIBMPI) -o $(TARGET)

clean:
	@echo "junk" > dummy.o
	@rm -f *.o  *~ $(TARGET_MPI) $(TARGET_SERIAL)
