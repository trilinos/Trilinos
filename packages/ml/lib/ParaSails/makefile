# Makefile for ParaSails
# Please modify for your local system.
# If you are using ESSL instead of LAPACK/BLAS, please define ESSL in Common.h

.SUFFIXES: .c .o

# this compiler automatically pulls in MPI
CC = /usr/local/mpich-1.2.0/bin/mpicc

C_COMPILE_FLAGS= -O -DPARASAILS_DEBUG -DPARASAILS_CG_PRINT
 
CFLAGS = \
 ${C_COMPILE_FLAGS}

RANLIB = ranlib

# need to link to LAPACK and BLAS
# linking fortran versions needs fortran libs
MPILIBFLAGS =
LIBFLAGS =
BLASLIBFLAGS = /home/casc/lib/lapack_solaris.a /home/casc/lib/blas_solaris.a
LDLIBFLAGS = -lF77 -lM77 -lsunmath -lm

LFLAGS =\
 -L. -lParaSails\
 ${BLASLIBFLAGS}\
 ${MPILIBFLAGS} ${LIBFLAGS} ${LDLIBFLAGS}

HEADERS =\
 Common.h\
 ConjGrad.h\
 DiagScale.h\
 Hash.h\
 LoadBal.h\
 Matrix.h\
 Mem.h\
 Numbering.h\
 OrderStat.h\
 ParaSails.h\
 PrunedRows.h\
 RowPatt.h\
 StoredRows.h

FILES =\
 ConjGrad.c\
 DiagScale.c\
 FlexGmres.c\
 Hash.c\
 LoadBal.c\
 Matrix.c\
 Mem.c\
 Numbering.c\
 OrderStat.c\
 ParaSails.c\
 PrunedRows.c\
 RowPatt.c\
 StoredRows.c

OBJS = ${FILES:.c=.o}

##################################################################
# Targets
##################################################################

all: libParaSails.a 

driver: driver.o libParaSails.a
	@echo  "Linking" $@ "... "
	${CC} -o $@ driver.o ${LFLAGS}

clean:
	@rm -f *.o libParaSails.a

##################################################################
# Rules
##################################################################

libParaSails.a: ${OBJS}
	@echo  "Building $@ ... "
	ar -rcu $@ ${OBJS}
	${RANLIB} $@

${OBJS}: ${HEADERS}

##################################################################
# Generic rules
##################################################################

.c.o:
	${CC} -o $@ -c ${CFLAGS} $<
