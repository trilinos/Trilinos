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

ML_C = ml_bdrypts.c      ml_elementagx.c   ml_get_basis.c \
       ml_grid.c         ml_gridfunc.c     ml_mapper.c \
       ml_gridagx.c      ml_pde.c          ml_setup.c

#=======================================================================
# ML include files
#=======================================================================

ML_INC = ml_bdrypts.h    ml_elementagx.h   ml_grid.h\
	 ml_gridfunc.h   ml_mapper.h       ml_gridagx.h \
         ml_pde.h 

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

ml_bdrypts.o : ml_bdrypts.c ml_bdrypts.h ../Include/ml_defs.h \
               ../Utils/ml_memory.h
	$(CC) -c $(CFLAGS) ml_bdrypts.c -o $@

ml_elementagx.o : ml_elementagx.c ml_elementagx.h ../Utils/ml_memory.h
	$(CC) -c $(CFLAGS) ml_elementagx.c -o $@

ml_get_basis.o : ml_get_basis.c ../Include/ml_defs.h ../Utils/ml_memory.h\
                 ml_gridfunc.h
	$(CC) -c $(CFLAGS) ml_get_basis.c -o $@

ml_grid.o : ml_grid.c ../Include/ml_defs.h ml_grid.h ml_gridfunc.h \
            ../Utils/ml_memory.h
	$(CC) -c $(CFLAGS) ml_grid.c -o $@

ml_gridagx.o : ml_gridagx.c ../Include/ml_defs.h ml_elementagx.h \
               ml_gridagx.h ../Utils/ml_intlist.h ../Utils/ml_memory.h
	$(CC) -c $(CFLAGS) ml_gridagx.c -o $@

ml_gridfunc.o : ml_gridfunc.c ../Include/ml_defs.h ml_gridfunc.h \
                ../Utils/ml_memory.h
	$(CC) -c $(CFLAGS) ml_gridfunc.c -o $@

ml_pde.o : ml_pde.c ../Comm/ml_comm.h ../Include/ml_defs.h \
           ../Utils/ml_memory.h ml_pde.h ../Main/mli_solver.h
	$(CC) -c $(CFLAGS) ml_pde.c -o $@

ml_mapper.o : ml_mapper.c ../Include/ml_defs.h ml_mapper.h \
              ../Utils/ml_memory.h
	$(CC) -c $(CFLAGS) ml_mapper.c -o $@

ml_setup.o : ml_setup.c ../Main/ml_1level.h ../FEGrid/ml_bdrypts.h \
             ../Comm/ml_comm.h ../Comm/ml_comminfoagx.h \
             ../Comm/ml_comminfoop.h ../Include/ml_defs.h \
             ../FEGrid/ml_grid.h ../FEGrid/ml_gridagx.h \
             ../FEGrid/ml_gridfunc.h ../Utils/ml_intlist.h \
             ../Utils/ml_memory.h ../Operator/ml_operator.h \
             ../Operator/ml_operatoragx.h ../Main/ml_struct.h \
             ../Utils/ml_utils.h ../FEGrid/ml_elementagx.h
	$(CC) -c $(CFLAGS) ml_setup.c -o $@

#=======================================================================

clean:
	@echo "junk" > dummy.o
	@rm -f *.o *.a core 

