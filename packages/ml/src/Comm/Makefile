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

INCLUDES= -I.. $(SUPERLU_INC_PATH) $(AZTEC_INC_PATH) $(MPI_INC_PATH) \
          -I../Include -I../Main -I../Utils -I../Krylov -I../Operator \
          -I../Comm -I../FEGrid -I../Smoother
LIB_PATHS=$(AZTEC_LIB_PATH) $(SUPERLU_LIB_PATH) $(MPI_LIB_PATH) 

CFLAGS=$(ARCH_FLAGS) $(DEFINES) $(INCLUDES) 

#=======================================================================
# ML source files
#=======================================================================

ML_C = ml_comm.c         ml_comminfoop.c   ml_comminfoagx.c \
       ml_exch_row.c

#=======================================================================
# ML include files
#=======================================================================

ML_INC = ml_comm.h       ml_comminfoop.h   ml_comminfoagx.h

#=======================================================================
# copy object files to library directory
#=======================================================================

ML_OBJ = $(ML_C:.c=.o)

tolib: $(ML_OBJ)
	@cp -f *.o ../lib

#=======================================================================
# dependencies for 'c' files
#=======================================================================
#include ../../etc/depends.ml

ml_comm.o : ml_comm.c ml_comm.h ../Include/ml_defs.h ../Utils/ml_memory.h
	$(CC) -c $(CFLAGS) ml_comm.c -o $@

ml_comminfoagx.o : ml_comminfoagx.c ml_comminfoagx.h ../Include/ml_defs.h\
                   ../Utils/ml_memory.h
	$(CC) -c $(CFLAGS) ml_comminfoagx.c -o $@

ml_comminfoop.o : ml_comminfoop.c ml_comm.h ../Include/ml_defs.h \
                  ../Utils/ml_memory.h
	$(CC) -c $(CFLAGS) ml_comminfoop.c -o $@

ml_exch_row.o : ml_exch_row.c ../Main/ml_1level.h \
                ml_comm.h ml_comminfoop.h ../Include/ml_defs.h \
                ../Operator/ml_operator.h ../Main/ml_struct.h \
                ../Utils/ml_utils.h ../Operator/ml_mat_formats.h
	$(CC) -c $(CFLAGS) ml_exch_row.c -o $@

#=======================================================================

clean:
	@echo "junk" > dummy.o
	@rm -f *.o *.a core 

