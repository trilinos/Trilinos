include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)
include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

DEFINES= -D$(TRILINOS_ARCH) $(EPETRA_ARCH_DEFINES) -D$(EPETRA_COMM)

KUNDERTSPARSE_INCLUDES = -I$(TRILINOS_HOME)/contrib/epetra_sparse
INCLUDES = $(ARCH_INCLUDES) $(EPETRA_INCLUDES) $(KUNDERTSPARSE_INCLUDES)

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES) $(INCLUDES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES) -DXyce_NEW_SPARSE
LDFLAGS=$(ARCH_LDFLAGS)

#=======================================================================
# EPETRA_SPARSE source files
#=======================================================================

LIBEPETRA_SPARSE = $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libepetra_sparse.a
LIBKUNDERT_SPARSE = $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)/libsparse.a

EPETRA_SPARSE_CC = Epetra_CrsKundertSparse.cpp
EPETRA_SPARSE_C = alloc.c   smplinkrows.c  spbuild.c   spoutput.c \
spsolve.c driver.c  spalloc.c	 spfactor.c  spsmp.c	 sputils.c


EPETRA_SPARSE_OBJ = $(EPETRA_SPARSE_CC:.cpp=.o) \
                    $(EPETRA_SPARSE_C:.c=.o)

$(LIBEPETRA_SPARSE): $(EPETRA_SPARSE_OBJ)
	mkdir -p $(TRILINOS_LIB_PATH)
	$(AR) $(ARFLAGS) $(LIBEPETRA_SPARSE) $(EPETRA_SPARSE_OBJ)
	$(RANLIB) $(LIBEPETRA_SPARSE)

clean:
	@rm -rf $(TRILINOS_TARGET)
	@rm -f *.o $(LIBEPETRA_SPARSE) *~ 
