include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)
include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

DEFINES= -D$(TRILINOS_ARCH) $(EPETRA_ARCH_DEFINES) -D$(EPETRA_COMM) $(EPETRA_DEFAULT_DEFINES)

ESI_INC = -I$(TRILINOS_HOME)/external/esi/cxx/include
EPETRA_INC = -I$(TRILINOS_HOME)/packages/epetra/src
AZTECOO_INC = -I$(TRILINOS_HOME)/packages/aztecoo/src

INCLUDES = $(ESI_INC) $(EPETRA_INC) $(AZTECOO_INC) $(ARCH_INCLUDES) \
           $(EPETRA_ESI_INCLUDES)

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) 
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
CCFLAGS=$(CXXFLAGS)
#=======================================================================
# EPETRA_ESI source files
#=======================================================================


EPETRA_ESI_SRC = \
	Epetra_ESI_inst.cpp

EPETRA_ESI_OBJ          = $(EPETRA_ESI_SRC:.cpp=.$(OBJ_EXT))

#$(LIBEPETRA_ESI): $(EPETRA_ESI_OBJ)
LIBEPETRA_ESI: $(EPETRA_ESI_OBJ)
	mkdir -p $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)
	$(AR) $(ARFLAGS) $(LIB_WIN_FLAGS)$(LIBEPETRA_ESI) $(EPETRA_ESI_OBJ)
# following ranlib commented out by rbl
#	$(RANLIB) $(LIBEPETRA_ESI)

clean:
	rm -rf $(TRILINOS_TARGET)
	rm -f *.$(OBJ_EXT)  $(LIBEPETRA_ESI) *~

