include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)
include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

DEFINES= -D$(TRILINOS_ARCH)

INCLUDES = $(ARCH_INCLUDES) $(TSF_INCLUDES)

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) 
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
CCFLAGS=$(CXXFLAGS)
#=======================================================================
# TSF source files
#=======================================================================


TSF_CC = TSF_Parameter.cc TSF_ParameterList.cc
TSF_C = 
TSF_F =

#=======================================================================
# TSF include files
#=======================================================================

TSF_INC =

TSF_OBJ          = $(TSF_CC:.cc=.o) $(TSF_C:.c=.o) $(TSF_F:.f=.o)

$(LIBTSF): $(TSF_OBJ)
	mkdir -p $(TRILINOS_HOME)/lib/$(TSF_TARGET)
	$(AR) $(ARFLAGS) $(LIBTSF) $(TSF_OBJ)
	$(RANLIB) $(LIBTSF)

#
# dependencies for 'f' files (none at this time)
#
#include ../../build/depends.trilinos

clean:
	@echo "junk" > dummy.o
	@rm -f *.o  $(LIBTSF) *~
