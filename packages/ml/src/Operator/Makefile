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

ML_C = ml_mat_formats.c     ml_matmat_mult.c     ml_op_utils.c\
       ml_operator.c        ml_operatoragx.c     ml_rap.c\
       ml_rap_utils.c

#=======================================================================
# ML include files
#=======================================================================

ML_INC = ml_mat_formats.h    ml_op_utils.h    ml_operator.h \
         ml_operatoragx.h    ml_rap.h

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

ml_mat_formats.o : ml_mat_formats.c ../Utils/ml_aztec_utils.h \
                   ../Main/ml_1level.h ../Comm/ml_comm.h\
                   ../Utils/ml_aztec_utils.h ../Comm/ml_comminfoagx.h \
                   ../Comm/ml_comminfoop.h ../FEGrid/ml_mapper.h\
                   ../Include/ml_defs.h ../Include/ml_include.h \
                   ml_mat_formats.h ../Utils/ml_memory.h ml_op_utils.h\
                   ml_operatoragx.h ml_rap.h ../Main/ml_struct.h \
                   ../Utils/ml_utils.h ../Utils/ml_vec.h ../FEGrid/ml_bdrypts.h
	@$(CC) -c $(CFLAGS) ml_mat_formats.c -o $@

ml_matmat_mult.o : ml_matmat_mult.c ../Main/ml_1level.h ../Comm/ml_comm.h\
                   ../Utils/ml_aztec_utils.h ../Comm/ml_comminfoagx.h \
                   ../Comm/ml_comminfoop.h ../FEGrid/ml_mapper.h\
                   ../Include/ml_defs.h ../Include/ml_include.h \
                   ml_mat_formats.h ../Utils/ml_memory.h ml_op_utils.h\
                   ml_operatoragx.h ml_rap.h ../Main/ml_struct.h \
                   ../Utils/ml_utils.h ../Utils/ml_vec.h ../FEGrid/ml_bdrypts.h
	@$(CC) -c $(CFLAGS) ml_matmat_mult.c -o $@

ml_op_utils.o : ml_op_utils.c ../Utils/ml_aztec_utils.h ../Comm/ml_comm.h\
                ../Comm/ml_comminfoagx.h ../Include/ml_defs.h \
                ../Include/ml_include.h ml_mat_formats.h ../Utils/ml_memory.h\
                ml_op_utils.h ml_operatoragx.h ml_rap.h \
                ../Main/ml_struct.h ../Utils/ml_utils.h ../Utils/ml_vec.h 
	@$(CC) -c $(CFLAGS) ml_op_utils.c -o $@

ml_operator.o : ml_operator.c ../Utils/ml_aztec_utils.h ../Comm/ml_comm.h\
                ../Comm/ml_comminfoagx.h ../Include/ml_defs.h \
                ../Main/ml_1level.h ../Include/ml_include.h ml_mat_formats.h \
                ../Utils/ml_memory.h ml_op_utils.h ml_operatoragx.h ml_rap.h \
                ../Main/ml_struct.h ../Utils/ml_utils.h ../Utils/ml_vec.h  \
                ../FEGrid/ml_bdrypts.h 
	@$(CC) -c $(CFLAGS) ml_operator.c -o $@

ml_operatoragx.o : ml_operatoragx.c ../Comm/ml_comm.h ml_operatoragx.h \
                   ../Comm/ml_comminfoagx.h ../Utils/ml_memory.h \
                   ../Main/ml_struct.h
	@$(CC) -c $(CFLAGS) ml_operatoragx.c -o $@

ml_rap.o : ml_rap.c ../Main/ml_1level.h ../FEGrid/ml_bdrypts.h \
           ../Comm/ml_comm.h ../Comm/ml_comminfoagx.h ../Comm/ml_comminfoop.h\
           ../Include/ml_defs.h ../FEGrid/ml_mapper.h ml_mat_formats.h \
           ml_operator.h ml_operatoragx.h ml_rap.h ../Main/ml_struct.h \
           ../Utils/ml_utils.h ../Utils/ml_vec.h
	@$(CC) -c $(CFLAGS) ml_rap.c -o $@

ml_rap_utils.o : ml_rap_utils.c ../Main/ml_1level.h ../FEGrid/ml_bdrypts.h \
                 ../Comm/ml_comm.h ../Comm/ml_comminfoop.h ../Include/ml_defs.h\
                 ../FEGrid/ml_mapper.h ml_mat_formats.h ml_operator.h \
                 ml_operatoragx.h ml_rap.h ../Main/ml_struct.h \
                 ../Utils/ml_utils.h ../Utils/ml_vec.h
	@$(CC) -c $(CFLAGS) ml_rap_utils.c -o $@

#=======================================================================

clean:
	@echo "junk" > dummy.o
	@rm -f *.o *.a core 

