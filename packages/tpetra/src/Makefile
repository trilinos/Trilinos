include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)
include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

DEFINES= -D$(TRILINOS_ARCH) $(TPETRA_ARCH_DEFINES)

INCLUDES = $(ARCH_INCLUDES) $(TPETRA_INCLUDES) 

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) 
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)
CCFLAGS=$(CXXFLAGS)
#=======================================================================
# TRILINOS_UTIL source files
#=======================================================================


TRILINOS_UTIL_CC = 

TRILINOS_UTIL_C = 
TRILINOS_UTIL_F =

#=======================================================================
# TRILINOS_UTIL include files
#=======================================================================

TRILINOS_UTIL_INC =

TRILINOS_UTIL_OBJ = $(TRILINOS_UTIL_CC:.cpp=.$(OBJ_EXT)) \
                    $(TRILINOS_UTIL_C:.c=.$(OBJ_EXT)) \
                    $(TRILINOS_UTIL_F:.for=.$(OBJ_EXT))

#$(LIBTRILINOS_UTIL): $(TRILINOS_UTIL_OBJ)
LIBTRILINOS_UTIL: $(TRILINOS_UTIL_OBJ)
	mkdir -p $(TRILINOS_HOME)/lib/$(TRILINOS_TARGET)
	$(AR) $(ARFLAGS) $(LIB_WIN_FLAGS)$(LIBTRILINOS_UTIL) $(TRILINOS_UTIL_OBJ)
	$(RANLIB) $(LIBTRILINOS_UTIL)

#
# dependencies for 'f' files (none at this time)
#
#include ../../build/depends.petra

clean:
	@echo "junk" > dummy.$(OBJ_EXT)
	@rm -f *.$(OBJ_EXT)  $(LIBTRILINOS_UTIL) *~
