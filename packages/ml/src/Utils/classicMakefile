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
          -I../Comm -I../FEGrid -I../Smoother -I../Coarsen
LIB_PATHS=$(AZTEC_LIB_PATH) $(SUPERLU_LIB_PATH) $(MPI_LIB_PATH) 

CFLAGS=$(ARCH_FLAGS) $(DEFINES) $(INCLUDES) 

#=======================================================================
# ML source files
#=======================================================================

ML_C = ml_aztec_utils.c   ml_intlist.c     ml_lapack.c \
       ml_memory.c        ml_rbm.c         ml_utils.c \
       ml_vec.c           ml_twogrid_analysis.c

#=======================================================================
# ML include files
#=======================================================================

ML_INC = ml_lapack.h     ml_aztec_utils.h    ml_intlist.h \
         ml_memory.h     ml_rbm.h            ml_utils.h \
         ml_vec.h        ml_twogrid_analysis.h

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

ml_aztec_utils.o : ml_aztec_utils.c ../Include/ml_defs.h ml_memory.h \
                   ../Main/ml_struct.h ../Operator/ml_operator.h \
                   ../Comm/ml_comminfoop.h
	$(CC) -c $(CFLAGS) ml_aztec_utils.c -o $@

ml_intlist.o : ml_intlist.c ../Include/ml_defs.h ml_intlist.h ml_memory.h
	$(CC) -c $(CFLAGS) ml_intlist.c -o $@

ml_lapack.o : ml_lapack.c ../Include/ml_defs.h
	$(CC) -c $(CFLAGS) ml_lapack.c -o $@

ml_memory.o : ml_memory.c ../Comm/ml_comm.h ../Include/ml_defs.h \
              ml_memory.h
	$(CC) -c $(CFLAGS) ml_memory.c -o $@

ml_rbm.o : ml_rbm.c ml_rbm.h 
	$(CC) -c $(CFLAGS) ml_rbm.c -o $@

ml_utils.o : ml_utils.c ../Comm/ml_comm.h ../Include/ml_defs.h ml_utils.h 
	$(CC) -c $(CFLAGS) ml_utils.c -o $@

ml_vec.o : ml_vec.c ../Comm/ml_comm.h ../Include/ml_defs.h ml_vec.h \
           ml_memory.h
	$(CC) -c $(CFLAGS) ml_vec.c -o $@

ml_twogrid_analysis.o : ml_twogrid_analysis.c ml_twogrid_analysis.h
	$(CC) -c $(CFLAGS) ml_twogrid_analysis.c -o $@

#=======================================================================

clean:
	@echo "junk" > dummy.o
	@rm -f *.o *.a core 

