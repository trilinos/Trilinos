# To run makefile:
#    1) set environment variable ARCH to sgi, sun, tflop, or pclinux.
#       Other machines require an appropriate makefile.$(ARCH) file.
#
#    2) MPI users:
#         a) set environment variable PACKAGES to "-DML_MPI"
#         b) edit makefile.$(ARCH) to use mpicc or 
#            set environment variables 
#                MPI_INC_PATH: -I/mpi_include_directories  
#                MPI_LIB_PATH: -L/mpi_lib_directories -lmpi
#
#    3) Aztec users:
#         a) set environment variable PACKAGES to "-DAZTEC"
#         b) set environment variables 
#             AZTEC_INC_PATH : -I/aztec_include directories
#             AZTEC_LIB_PATH : -L/aztec_lib_directories -laztec
#
#    4) SuperLU users:
#         a) set environment variable PACKAGES to "-DSUPERLU"
#         b) set environment variables 
#             SUPERLU_INC_PATH : -I/superlu_include directories
#             SUPERLU_LIB_PATH : -L/superlu_lib_directories -lsuperlu
#
include ../../etc/makefile.$(ARCH)
DEFINES= -D$(ARCH) $(PACKAGES) 

INCLUDES=-I../src $(SUPERLU_INC_PATH) $(AZTEC_INC_PATH) $(MPI_INC_PATH)\
         $(PARASAILS_INC_PATH) \
         -I../Coarsen -I../Comm -I../FEGrid -I../Operator -I../Main \
         -I../Utils -I../Krylov -I../Smoother -I../Include -I..

LIB_PATHS=$(AZTEC_LIB_PATH) $(SUPERLU_LIB_PATH) $(PARASAILS_LIB_PATH) $(MPI_LIB_PATH) 

CFLAGS=$(ARCH_FLAGS) $(DEFINES) $(INCLUDES) 

#=======================================================================

#=======================================================================

APPLICATION_C = ml_example2d.c
APPLICATION_OBJ = $(APPLICATION_C:.c=.o)
APPLICATION     = $(APPLICATION_C:.c= )

$(APPLICATION): $(APPLICATION_OBJ) ../lib/libml.a
	(cd ..; make)
	echo "Linking"; \
	$(LINKER) $(APPLICATION_OBJ) $(ARCH_FLAGS) $(DEFINES) \
	   $(INCLUDES) -o $(APPLICATION) ../lib/libml.a $(LIB_PATHS) \
	   $(ARCH_LIBS) -lm

driver: driver.o ../lib/libml.a
	(cd ..; make)
	echo "Linking"; \
	$(LINKER) driver.o $(ARCH_FLAGS) $(DEFINES) \
           $(INCLUDES) -o driver ../lib/libml.a $(LIB_PATHS) \
           $(ARCH_LIBS) -lm

#
# dependencies for 'c' files
#
#include ../etc/depends.ml

ml_example2d.o : ../../examples/ml_example2d.c ../Coarsen/ml_aggregate.h \
        ../Utils/ml_aztec_utils.h ../Comm/ml_comm.h ../Include/ml_defs.h\
        ../Include/ml_include.h ../Utils/ml_memory.h ../Include/ml_include.h\
        ../Main/ml_struct.h ../Utils/ml_utils.h
	$(CC) -c $(CFLAGS) ../../examples/ml_example2d.c -o $@

driver.o : ../../examples/ml_example2d.c ../Coarsen/ml_aggregate.h \
        ../Comm/ml_comm.h ../Include/ml_defs.h\
        ../Include/ml_include.h ../Utils/ml_memory.h ../Include/ml_include.h\
        ../Main/ml_struct.h ../Utils/ml_utils.h
	$(CC) -c $(CFLAGS) ../Main/driver.c -o $@

clean:
	@echo "junk" > dummy.o
	@rm -f *.o $(APPLICATION) core PI*

