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

ML_C = ml_krylov.c   ml_cg.c    ml_gmres.c   ml_bicgstabl.c

#=======================================================================
# ML include files
#=======================================================================

ML_INC = ml_krylov.h   ml_cg.h    ml_gmres.h   ml_bicgstabl.h

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

ml_krylov.o : ml_krylov.c ../Comm/ml_comm.h ml_krylov.h \
              ../Operator/ml_operator.h ../Utils/ml_memory.h \
              ../FEGrid/ml_bdrypts.h ../Main/ml_1level.h \
              ../Utils/ml_vec.h ../Include/ml_defs.h
	$(CC) -c $(CFLAGS) ml_krylov.c -o $@

ml_cg.o : ml_cg.c ../Comm/ml_comm.h ml_cg.h ml_krylov.h \
          ../Operator/ml_operator.h ../Utils/ml_memory.h \
          ../FEGrid/ml_bdrypts.h ../Main/ml_1level.h ../Utils/ml_vec.h\
          ../Include/ml_defs.h
	$(CC) -c $(CFLAGS) ml_cg.c -o $@

ml_gmres.o : ml_gmres.c ml_gmres.h ml_krylov.h ../Comm/ml_comm.h \
             ../Operator/ml_operator.h ../Utils/ml_memory.h \
             ../FEGrid/ml_bdrypts.h ../Main/ml_1level.h ../Utils/ml_vec.h\
             ../Include/ml_defs.h
	$(CC) -c $(CFLAGS) ml_gmres.c -o $@

ml_bicgstabl.o : ml_bicgstabl.c ml_bicgstabl.h ml_krylov.h \
                 ../Comm/ml_comm.h ../Operator/ml_operator.h \
                 ../Utils/ml_memory.h ../FEGrid/ml_bdrypts.h \
                 ../Main/ml_1level.h ../Utils/ml_vec.h ../Include/ml_defs.h
	$(CC) -c $(CFLAGS) ml_bicgstabl.c -o $@

#=======================================================================

clean:
	@echo "junk" > dummy.o
	@rm -f *.o *.a core 

