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

ML_C = ml_struct.c  mli_solver.c  ml_seg_precond.c 

#=======================================================================
# ML include files
#=======================================================================

ML_INC = ml_struct.h   mli_solver.h   ml_seg_precond.h

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

ml_struct.o : ml_struct.c ml_1level.h ../Coarsen/ml_aggregate.h\
              ../FEGrid/ml_bdrypts.h ../Coarsen/ml_check.h ../Comm/ml_comm.h\
              ../Comm/ml_comminfoagx.h ../Comm/ml_comminfoop.h \
              ../Smoother/ml_csolve.h ../Include/ml_defs.h \
              ../FEGrid/ml_grid.h ../FEGrid/ml_gridfunc.h\
              ../Utils/ml_lapack.h ../Coarsen/ml_agg_genP.h \
              ../Krylov/ml_krylov.h ../FEGrid/ml_mapper.h \
              ../Operator/ml_mat_formats.h ../Utils/ml_memory.h \
              ../Operator/ml_operator.h ../Operator/ml_operatoragx.h \
              ../Operator/ml_rap.h ../Smoother/ml_smoother.h \
              ../Smoother/ml_solver.h ml_struct.h ../Utils/ml_utils.h\
              ../Utils/ml_vec.h
	@$(CC) -c $(CFLAGS) ml_struct.c -o $@

mli_solver.o : mli_solver.c ml_1level.h ../Coarsen/ml_aggregate.h\
              ../FEGrid/ml_bdrypts.h ../Coarsen/ml_check.h ../Comm/ml_comm.h\
              ../Comm/ml_comminfoagx.h ../Comm/ml_comminfoop.h \
              ../Smoother/ml_csolve.h ../Include/ml_defs.h \
              ../FEGrid/ml_grid.h ../FEGrid/ml_gridfunc.h\
              ../Utils/ml_lapack.h ../Coarsen/ml_agg_genP.h \
              ../Krylov/ml_krylov.h ../FEGrid/ml_mapper.h \
              ../Operator/ml_mat_formats.h ../Utils/ml_memory.h \
              ../Operator/ml_operator.h ../Operator/ml_operatoragx.h \
              ../Operator/ml_rap.h ../Smoother/ml_smoother.h \
              ../Smoother/ml_solver.h ml_struct.h ../Utils/ml_utils.h\
              ../Utils/ml_vec.h mli_solver.h
	@$(CC) -c $(CFLAGS) mli_solver.c -o $@

ml_seg_precond.o : ml_seg_precond.c ml_1level.h ../Coarsen/ml_aggregate.h\
              ../FEGrid/ml_bdrypts.h ../Coarsen/ml_check.h ../Comm/ml_comm.h\
              ../Comm/ml_comminfoagx.h ../Comm/ml_comminfoop.h \
              ../Smoother/ml_csolve.h ../Include/ml_defs.h \
              ../FEGrid/ml_grid.h ../FEGrid/ml_gridfunc.h\
              ../Utils/ml_lapack.h ../Coarsen/ml_agg_genP.h \
              ../Krylov/ml_krylov.h ../FEGrid/ml_mapper.h \
              ../Operator/ml_mat_formats.h ../Utils/ml_memory.h \
              ../Operator/ml_operator.h ../Operator/ml_operatoragx.h \
              ../Operator/ml_rap.h ../Smoother/ml_smoother.h \
              ../Smoother/ml_solver.h ml_struct.h ../Utils/ml_utils.h\
              ../Utils/ml_vec.h ml_seg_precond.h
	@$(CC) -c $(CFLAGS) ml_seg_precond.c -o $@

#=======================================================================

clean:
	@echo "junk" > dummy.o
	@rm -f *.o *.a core 

