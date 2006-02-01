#-------------------------------------------------------------------------------
# Makefile for paraklete
#
# PARAKLETE version 0.1: parallel sparse LU factorization.  May 13, 2005
# Copyright (C) 2005, Univ. of Florida.  Author: Timothy A. Davis
# See License.txt for the Version 2.1 of the GNU Lesser General Public License
# http://www.cise.ufl.edu/research/sparse
# 
#-------------------------------------------------------------------------------

all: pk
# all: parasym.mexglx  parafac.mexglx

# MATLAB:
#CC = mex -g

# valgrind, no MPI, pretend to have 5 processes, debug
# CC = gcc -Wall -W -Wshadow -Wmissing-prototypes -Wstrict-prototypes \
 	-Wredundant-decls -Wnested-externs -Wdisabled-optimization \
 	-pedantic -ansi -O -g -DNMPI -DHACK=5

# valgrind, no MPI, one process, no debug
# CC = gcc -Wall -W -Wshadow -Wmissing-prototypes -Wstrict-prototypes \
 	-Wredundant-decls -Wnested-externs -Wdisabled-optimization \
 	-pedantic -ansi -O3 -DNMPI

# MPI:
# CC = mpicc -Wall -W -Wshadow -Wmissing-prototypes -Wstrict-prototypes \
 	-Wredundant-decls -Wnested-externs -Wdisabled-optimization \
 	-pedantic -ansi -O -g

# MPI, no debug
CC = mpicc -Wall -W -Wshadow -Wmissing-prototypes -Wstrict-prototypes \
 	-Wredundant-decls -Wnested-externs -Wdisabled-optimization \
 	-pedantic -ansi -O3

# MATLAB:
#LIB = -lmetis -lsecond
# plain blas:
# LIB = -lmetis -lsecond -lm -llapack_plain -lblas_plain -lg2c
# LIB = -lmetis -lsecond -lm -llapack -lgoto -lg2c -lefence -lpthread
LIB = -L. -lmetis -lsecond -lm -llapack -lgoto -lxerbla -lg2c
LIB = -L. -lmetis

PSRC = paraklete_analyze.c  \
	paraklete_factorize.c \
	paraklete_kernel.c \
	paraklete_node.c \
	paraklete_solve.c \
	paraklete_lsolve_node.c \
	paraklete_usolve_node.c

PINC = paraklete.h

POBJ = paraklete_analyze.o \
	paraklete_factorize.o \
	paraklete_kernel.o \
	paraklete_node.o \
	paraklete_solve.o \
	paraklete_lsolve_node.o \
	paraklete_usolve_node.o

CMINC = ../Cholmod/Include/cholmod_check.h \
	../Cholmod/Include/cholmod_core.h \
	../Cholmod/Include/cholmod.h \
	../Cholmod/Include/cholmod_internal.h \
	../Cholmod/Include/cholmod_partition.h

AMDSRC = ../AMD/Source/amd_1.c \
	../AMD/Source/amd_2.c \
	../AMD/Source/amd_aat.c \
	../AMD/Source/amd_control.c \
	../AMD/Source/amd_defaults.c \
	../AMD/Source/amd_dump.c \
	../AMD/Source/amd_info.c \
	../AMD/Source/amd_order.c \
	../AMD/Source/amd_postorder.c \
	../AMD/Source/amd_post_tree.c \
	../AMD/Source/amd_preprocess.c \
	../AMD/Source/amd_valid.c \
	../AMD/Include/amd.h \
	../AMD/Source/amd_internal.h

AMDOBJ = amd_1.o \
	amd_2.o \
	amd_aat.o \
	amd_control.o \
	amd_defaults.o \
	amd_dump.o \
	amd_info.o \
	amd_order.o \
	amd_postorder.o \
	amd_post_tree.o \
	amd_preprocess.o \
	amd_valid.o

COLAMDSRC = ../colamd/colamd.c ../colamd/colamd.h

COLAMDOBJ = colamd.o

CCOLAMDSRC = ../CColamd/ccolamd.c \
	../CColamd/ccolamd_apply_order.c \
	../CColamd/ccolamd_fsize.c \
	../CColamd/ccolamd_post_tree.c \
	../CColamd/ccolamd_postorder.c \
	../CColamd/ccolamd.h \
	../CColamd/ccolamd_internal.h

CCOLAMDOBJ = ccolamd.o \
	ccolamd_apply_order.o \
	ccolamd_fsize.o \
	ccolamd_post_tree.o \
	ccolamd_postorder.o

SRC = \
	../Cholmod/Core/cholmod_common.c \
	../Cholmod/Core/cholmod_dense.c \
	../Cholmod/Core/cholmod_factor.c \
	../Cholmod/Core/cholmod_memory.c \
	../Cholmod/Core/cholmod_postorder.c \
	../Cholmod/Core/cholmod_sparse.c \
	../Cholmod/Core/cholmod_band.c \
	../Cholmod/Core/cholmod_copy.c \
	../Cholmod/Core/cholmod_triplet.c \
	../Cholmod/Core/cholmod_error.c \
	../Cholmod/Core/cholmod_aat.c \
	../Cholmod/Core/cholmod_add.c \
	../Cholmod/Check/cholmod_check.c \
	../Cholmod/Cholesky/cholmod_amd.c \
	../Cholmod/Cholesky/cholmod_analyze.c \
	../Cholmod/Cholesky/cholmod_colamd.c \
	../Cholmod/Cholesky/cholmod_etree.c \
	../Cholmod/Cholesky/cholmod_factorize.c \
	../Cholmod/Cholesky/cholmod_solve.c \
	../Cholmod/Cholesky/cholmod_lsolve.c \
	../Cholmod/Cholesky/cholmod_ltsolve.c \
	../Cholmod/Cholesky/cholmod_resymbol.c \
	../Cholmod/Cholesky/cholmod_rowcolcounts.c \
	../Cholmod/Cholesky/cholmod_rowfac.c \
	../Cholmod/Partition/cholmod_ccolamd.c \
	../Cholmod/Partition/cholmod_csymamd.c \
	../Cholmod/Partition/cholmod_metis.c \
	../Cholmod/Partition/cholmod_nd.c

CMOBJ = \
	cholmod_common.o \
	cholmod_dense.o \
	cholmod_factor.o \
	cholmod_memory.o \
	cholmod_postorder.o \
	cholmod_sparse.o \
	cholmod_band.o \
	cholmod_copy.o \
	cholmod_triplet.o \
	cholmod_error.o \
	cholmod_aat.o \
	cholmod_add.o \
	cholmod_check.o \
	cholmod_amd.o \
	cholmod_analyze.o \
	cholmod_colamd.o \
	cholmod_etree.o \
	cholmod_factorize.o \
	cholmod_resymbol.o \
	cholmod_rowcolcounts.o \
	cholmod_rowfac.o \
	cholmod_solve.o \
	cholmod_ccolamd.o \
	cholmod_csymamd.o \
	cholmod_metis.o \
	cholmod_nd.o

# These files are #include'd in other files, so there are no *.o files:
UPSRC = ../Cholmod/Modify/cholmod_updown_numeric.c ../Cholmod/Modify/cholmod_updown_numkr.c
SOLSRC = ../Cholmod/Cholesky/cholmod_lsolve.c ../Cholmod/Cholesky/cholmod_ltsolve.c

$(AMDOBJ): $(AMDSRC)

$(CCOLAMDOBJ): $(CCOLAMDSRC)

$(CMOBJ): $(CMINC)

$(POBJ): $(PSRC) $(PINC) Makefile


I = -I../AMD/Include -I../AMD/Source -I../colamd \
	-I../metis-4.0/Lib -I../CColamd -I../Cholmod/Include

ALLOBJ = $(AMDOBJ) $(CCOLAMDOBJ) $(CMOBJ) $(COLAMDOBJ) $(POBJ)

#-------------------------------------------------------------------------------

pk: pk.c $(ALLOBJ) memory.c Makefile pk.h
	$(CC) $(I) -o pk pk.c memory.c $(ALLOBJ) $(LIB)

parasym.mexglx: parasym_mex.c $(ALLOBJ)
	$(CC) $(I) -o parasym parasym_mex.c $(ALLOBJ) $(LIB)

parafac.mexglx: parafac_mex.c $(ALLOBJ)
	$(CC) $(I) -o parafac parafac_mex.c $(ALLOBJ) $(LIB)

.c.o:
	$(CC) -c $(I) $*.c

#-------------------------------------------------------------------------------

amd_1.o: ../AMD/Source/amd_1.c
	$(CC) -c $(I) $<

amd_2.o: ../AMD/Source/amd_2.c
	$(CC) -c $(I) $<

amd_aat.o: ../AMD/Source/amd_aat.c
	$(CC) -c $(I) $<

amd_control.o: ../AMD/Source/amd_control.c
	$(CC) -c $(I) $<

amd_defaults.o: ../AMD/Source/amd_defaults.c
	$(CC) -c $(I) $<

amd_dump.o: ../AMD/Source/amd_dump.c
	$(CC) -c $(I) $<

amd_info.o: ../AMD/Source/amd_info.c
	$(CC) -c $(I) $<

amd_order.o: ../AMD/Source/amd_order.c
	$(CC) -c $(I) $<

amd_postorder.o: ../AMD/Source/amd_postorder.c
	$(CC) -c $(I) $<

amd_post_tree.o: ../AMD/Source/amd_post_tree.c
	$(CC) -c $(I) $<

amd_preprocess.o: ../AMD/Source/amd_preprocess.c
	$(CC) -c $(I) $<

amd_valid.o: ../AMD/Source/amd_valid.c
	$(CC) -c $(I) $<

#-------------------------------------------------------------------------------

colamd.o: ../colamd/colamd.c
	$(CC) -c $(I) $<

#-------------------------------------------------------------------------------

ccolamd.o: ../CColamd/ccolamd.c
	$(CC) -c $(I) $<

ccolamd_apply_order.o: ../CColamd/ccolamd_apply_order.c
	$(CC) -c $(I) $<

ccolamd_fsize.o: ../CColamd/ccolamd_fsize.c
	$(CC) -c $(I) $<

ccolamd_post_tree.o: ../CColamd/ccolamd_post_tree.c
	$(CC) -c $(I) $<

ccolamd_postorder.o: ../CColamd/ccolamd_postorder.c
	$(CC) -c $(I) $<

#-------------------------------------------------------------------------------

cholmod_check.o: ../Cholmod/Check/cholmod_check.c
	$(CC) -c $(I) $<

#-------------------------------------------------------------------------------

cholmod_common.o: ../Cholmod/Core/cholmod_common.c
	$(CC) -c $(I) $<

cholmod_dense.o: ../Cholmod/Core/cholmod_dense.c
	$(CC) -c $(I) $<

cholmod_factor.o: ../Cholmod/Core/cholmod_factor.c
	$(CC) -c $(I) $<

cholmod_memory.o: ../Cholmod/Core/cholmod_memory.c
	$(CC) -c $(I) $<

cholmod_postorder.o: ../Cholmod/Core/cholmod_postorder.c
	$(CC) -c $(I) $<

cholmod_sparse.o: ../Cholmod/Core/cholmod_sparse.c
	$(CC) -c $(I) $<

cholmod_band.o: ../Cholmod/Core/cholmod_band.c
	$(CC) -c $(I) $<

cholmod_copy.o: ../Cholmod/Core/cholmod_copy.c
	$(CC) -c $(I) $<

cholmod_triplet.o: ../Cholmod/Core/cholmod_triplet.c
	$(CC) -c $(I) $<

cholmod_error.o: ../Cholmod/Core/cholmod_error.c
	$(CC) -c $(I) $<

cholmod_aat.o: ../Cholmod/Core/cholmod_aat.c
	$(CC) -c $(I) $<

cholmod_add.o: ../Cholmod/Core/cholmod_add.c
	$(CC) -c $(I) $<

#-------------------------------------------------------------------------------

cholmod_amd.o: ../Cholmod/Cholesky/cholmod_amd.c
	$(CC) -c $(I) $<

cholmod_analyze.o: ../Cholmod/Cholesky/cholmod_analyze.c
	$(CC) -c $(I) $<

cholmod_colamd.o: ../Cholmod/Cholesky/cholmod_colamd.c
	$(CC) -c $(I) $<

cholmod_etree.o: ../Cholmod/Cholesky/cholmod_etree.c
	$(CC) -c $(I) $<

cholmod_factorize.o: ../Cholmod/Cholesky/cholmod_factorize.c
	$(CC) -c $(I) $<

cholmod_resymbol.o: ../Cholmod/Cholesky/cholmod_resymbol.c
	$(CC) -c $(I) $<

cholmod_rowcolcounts.o: ../Cholmod/Cholesky/cholmod_rowcolcounts.c
	$(CC) -c $(I) $<

cholmod_solve.o: ../Cholmod/Cholesky/cholmod_solve.c $(SOLSRC)
	$(CC) -c $(I) $<

cholmod_rowfac.o: ../Cholmod/Cholesky/cholmod_rowfac.c
	$(CC) -c $(I) $<

#------------------------------------------------------------------------------- 
cholmod_ccolamd.o: ../Cholmod/Partition/cholmod_ccolamd.c
	$(CC) -c $(I) $<

cholmod_csymamd.o: ../Cholmod/Partition/cholmod_csymamd.c
	$(CC) -c $(I) $<

cholmod_metis.o: ../Cholmod/Partition/cholmod_metis.c
	$(CC) -c $(I) $<

cholmod_nd.o: ../Cholmod/Partition/cholmod_nd.c
	$(CC) -c $(I) $<

#-------------------------------------------------------------------------------

cholmod_horzcat.o: ../Cholmod/MatrixOps/cholmod_horzcat.c
	$(CC) -c $(I) $<

cholmod_norm.o: ../Cholmod/MatrixOps/cholmod_norm.c
	$(CC) -c $(I) $<

cholmod_scale.o: ../Cholmod/MatrixOps/cholmod_scale.c
	$(CC) -c $(I) $<

cholmod_drop.o: ../Cholmod/MatrixOps/cholmod_drop.c
	$(CC) -c $(I) $<

cholmod_sdmult.o: ../Cholmod/MatrixOps/cholmod_sdmult.c
	$(CC) -c $(I) $<

cholmod_ssmult.o: ../Cholmod/MatrixOps/cholmod_ssmult.c
	$(CC) -c $(I) $<

cholmod_submatrix.o: ../Cholmod/MatrixOps/cholmod_submatrix.c
	$(CC) -c $(I) $<

cholmod_vertcat.o: ../Cholmod/MatrixOps/cholmod_vertcat.c
	$(CC) -c $(I) $<

#-------------------------------------------------------------------------------

cholmod_rowadd.o: ../Cholmod/Modify/cholmod_rowadd.c
	$(CC) -c $(I) $<

cholmod_rowdel.o: ../Cholmod/Modify/cholmod_rowdel.c
	$(CC) -c $(I) $<

cholmod_updown.o: ../Cholmod/Modify/cholmod_updown.c $(UPSRC)
	$(CC) -c $(I) $<

#-------------------------------------------------------------------------------

cholmod_super_numeric.o: ../Cholmod/Supernodal/cholmod_super_numeric.c
	$(CC) -c $(I) $<

cholmod_super_symbolic.o: ../Cholmod/Supernodal/cholmod_super_symbolic.c
	$(CC) -c $(I) $<

cholmod_super_solve.o: ../Cholmod/Supernodal/cholmod_super_solve.c
	$(CC) -c $(I) $<

#-------------------------------------------------------------------------------

clean:
	'rm' -f *.o

purge: clean
	'rm' *.mex* pk
