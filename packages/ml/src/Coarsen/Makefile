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

ML_C = ml_ggraph.c       ml_aggregate.c    ml_agg_genP.c \
       ml_agg_MIS.c      ml_agg_coupled.c  ml_agg_uncoupled.c \
       ml_agg_dd.c       ml_check.c        ml_amg_genP.c \
       ml_amg.c          ml_amg_MIS.c

#=======================================================================
# ML include files
#=======================================================================

ML_INC = ml_ggraph.h         ml_aggregate.h      ml_agg_genP.h\
	 ml_check.h          ml_amg.h

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

ml_aggregate.o : ml_aggregate.c ml_aggregate.h ../Comm/ml_comm.h \
             ../Include/ml_defs.h ../Utils/ml_lapack.h ../Utils/ml_memory.h\
             ../Operator/ml_operator.h
	@$(CC) -c $(CFLAGS) ml_aggregate.c -o $@

ml_check.o : ml_check.c ../Main/ml_1level.h ../FEGrid/ml_bdrypts.h \
             ml_check.h ../Comm/ml_comm.h ../Comm/ml_comminfoagx.h \
             ../Comm/ml_comminfoop.h ../Smoother/ml_csolve.h \
             ../Include/ml_defs.h ../FEGrid/ml_grid.h ../FEGrid/ml_gridfunc.h\
             ../Krylov/ml_krylov.h ../FEGrid/ml_mapper.h \
             ../Operator/ml_mat_formats.h ../Utils/ml_memory.h \
             ../Operator/ml_operator.h ../Operator/ml_operatoragx.h \
             ../Operator/ml_rap.h ../Smoother/ml_smoother.h \
             ../Smoother/ml_solver.h ../Main/ml_struct.h ../Utils/ml_utils.h\
             ../Utils/ml_vec.h
	@$(CC) -c $(CFLAGS) ml_check.c -o $@

ml_agg_genP.o : ml_agg_genP.c ml_aggregate.h \
        ../Utils/ml_aztec_utils.h ../Krylov/ml_bicgstabl.h ../Krylov/ml_cg.h\
        ../Comm/ml_comm.h ../Comm/ml_comminfoagx.h ../Include/ml_defs.h \
        ../FEGrid/ml_elementagx.h ml_ggraph.h ../Krylov/ml_gmres.h \
        ../FEGrid/ml_gridagx.h ../FEGrid/ml_gridfunc.h ../Include/ml_include.h\
        ../Utils/ml_intlist.h ../Krylov/ml_krylov.h ../Operator/ml_mat_formats.h\
        ../Utils/ml_memory.h ../Operator/ml_op_utils.h ../Operator/ml_operator.h\
        ../Operator/ml_operatoragx.h ../Operator/ml_rap.h ../Utils/ml_rbm.h\
        ml_agg_genP.h ../Smoother/ml_solver.h ../Main/ml_struct.h\
        ../Utils/ml_utils.h ../Utils/ml_vec.h ../Main/mli_solver.h
	@$(CC) -c $(CFLAGS) ml_agg_genP.c -o $@

ml_ggraph.o : ml_ggraph.c ../Comm/ml_comm.h ../Comm/ml_comminfoop.h \
        ../Include/ml_defs.h ml_ggraph.h ../FEGrid/ml_gridfunc.h \
        ../Operator/ml_mat_formats.h ../Utils/ml_memory.h \
        ../Operator/ml_operator.h
	@$(CC) -c $(CFLAGS) ml_ggraph.c -o $@

ml_agg_uncoupled.o : ml_agg_uncoupled.c ml_aggregate.h ../Comm/ml_comm.h \
        ../Comm/ml_comminfoop.h ../Include/ml_defs.h ../Operator/ml_operator.h\
        ../Operator/ml_mat_formats.h ../Utils/ml_memory.h
	@$(CC) -c $(CFLAGS) ml_agg_uncoupled.c -o $@

ml_agg_coupled.o : ml_agg_coupled.c ml_aggregate.h ../Comm/ml_comm.h \
        ../Comm/ml_comminfoop.h ../Include/ml_defs.h ../Operator/ml_operator.h\
        ../Operator/ml_mat_formats.h ../Utils/ml_memory.h
	@$(CC) -c $(CFLAGS) ml_agg_coupled.c -o $@

ml_agg_MIS.o : ml_agg_coupled.c ml_aggregate.h ml_ggraph.h ../Comm/ml_comm.h \
        ../Comm/ml_comminfoop.h ../Include/ml_defs.h ../Operator/ml_operator.h\
        ../Operator/ml_mat_formats.h ../Utils/ml_memory.h
	@$(CC) -c $(CFLAGS) ml_agg_MIS.c -o $@

ml_agg_dd.o : ml_agg_dd.c ml_aggregate.h ../Comm/ml_comm.h \
        ../Comm/ml_comminfoop.h ../Include/ml_defs.h ../Operator/ml_operator.h\
        ../Operator/ml_mat_formats.h ../Utils/ml_memory.h
	@$(CC) -c $(CFLAGS) ml_agg_dd.c -o $@

ml_amg.o : ml_amg.c ml_amg.h ../Comm/ml_comm.h \
             ../Include/ml_defs.h ../Utils/ml_lapack.h ../Utils/ml_memory.h\
             ../Operator/ml_operator.h
	@$(CC) -c $(CFLAGS) ml_amg.c -o $@

ml_amg_MIS.o : ml_amg_MIS.c ml_amg.h ../Comm/ml_comm.h \
             ../Include/ml_defs.h ../Utils/ml_lapack.h ../Utils/ml_memory.h\
             ../Operator/ml_operator.h
	@$(CC) -c $(CFLAGS) ml_amg_MIS.c -o $@

ml_amg_genP.o : ml_amg_genP.c ml_amg.h ../Comm/ml_comm.h \
             ../Include/ml_defs.h ../Utils/ml_lapack.h ../Utils/ml_memory.h\
             ../Operator/ml_operator.h
	@$(CC) -c $(CFLAGS) ml_amg_genP.c -o $@

#=======================================================================

clean:
	@echo "junk" > dummy.o
	@rm -f *.o *.a core 

