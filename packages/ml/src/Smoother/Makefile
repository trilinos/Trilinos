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
          $(PARASAILS_INC_PATH) \
          -I../Include -I../Main -I../Utils -I../Krylov -I../Operator \
          -I../Comm -I../FEGrid -I../Smoother -I../Coarsen
LIB_PATHS=$(AZTEC_LIB_PATH) $(SUPERLU_LIB_PATH) $(MPI_LIB_PATH) 

CFLAGS=$(ARCH_FLAGS) $(DEFINES) $(INCLUDES) 

#=======================================================================
# ML source files
#=======================================================================

ML_C = ml_csolve.c     ml_smoother.c     ml_solver.c\
       ml_superlu.c    ml_xxt.c

#=======================================================================
# ML include files
#=======================================================================

ML_INC = ml_csolve.h    ml_solver.h    ml_smoother.h

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

ml_csolve.o : ml_csolve.c ../Utils/ml_aztec_utils.h ../Comm/ml_comm.h\
              ../Comm/ml_comminfoagx.h ../Include/ml_defs.h \
              ../Operator/ml_mat_formats.h \
              ../Utils/ml_memory.h ../Operator/ml_op_utils.h \
              ../Main/ml_struct.h ../Utils/ml_utils.h ml_csolve.h
	$(CC) -c $(CFLAGS) ml_csolve.c -o $@

ml_solver.o : ml_solver.c ../Utils/ml_aztec_utils.h ../Comm/ml_comm.h\
              ../Comm/ml_comminfoagx.h ../Include/ml_defs.h \
              ../Operator/ml_mat_formats.h \
              ../Utils/ml_memory.h ../Operator/ml_op_utils.h \
              ../Main/ml_struct.h ../Utils/ml_utils.h ml_solver.h
	$(CC) -c $(CFLAGS) ml_solver.c -o $@

ml_smoother.o : ml_smoother.c ../Main/ml_1level.h ../Include/ml_defs.h \
                ../Utils/ml_lapack.h ../Utils/ml_memory.h ml_smoother.h \
                ../Operator/ml_operator.h
	$(CC) -c $(CFLAGS) ml_smoother.c -o $@

ml_superlu.o : ml_superlu.c ../Comm/ml_comm.h ml_solver.h \
               ../Include/ml_defs.h ../Operator/ml_mat_formats.h \
               ../Utils/ml_memory.h ../Main/ml_struct.h ../Utils/ml_utils.h 
	$(CC) -c $(CFLAGS) ml_superlu.c -o $@

ml_xxt.o : ml_xxt.c ../Comm/ml_comm.h ../Comm/ml_comminfoop.h \
           ../Include/ml_defs.h \
           ../Operator/ml_mat_formats.h ../Utils/ml_memory.h \
           ../Operator/ml_op_utils.h ../Main/ml_struct.h ../Utils/ml_utils.h
	$(CC) -c $(CFLAGS) ml_xxt.c -o $@


#=======================================================================

clean:
	@echo "junk" > dummy.o
	@rm -f *.o *.a core 

