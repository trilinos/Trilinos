# To run makefile:
#    1) set environment variable TRILINOS_ARCH to sgi, sun, tflop, or 
#    pclinux.
#       Other machines require an appropriate 
#       makefile.$(TRILINOS_ARCH) file.
#    2) Set TRILINOS_COMM to SERIAL or MPI
#    3) (Optional) Set TRILINOS_ID to make unique version for same
#       architecture and communication mode.
#
#    4) Make the archive $(LIBY12M) by typing 'make'.
#

include $(TRILINOS_HOME)/build/TRILINOS_TARGET_DEFS
TRILINOS_TARGET = $(TRILINOS_ARCH).$(TRILINOS_COMM)$(TRILINOS_ID)

include $(TRILINOS_HOME)/build/makefile.$(TRILINOS_TARGET)

DEFINES= -D$(TRILINOS_ARCH)

INCLUDES =
LIB_PATHS=

CFLAGS=$(ARCH_CFLAGS) $(DEFINES) $(INCLUDES)
FFLAGS=$(ARCH_FFLAGS) $(DEFINES)
CXXFLAGS=$(ARCH_CXXFLAGS) $(DEFINES) $(INCLUDES)

#=======================================================================
# Y12M source files
#=======================================================================

Y12M_F = y12mae.f  y12maf.f  y12mfe.f

#=======================================================================
# Y12M include files
#=======================================================================

Y12M_INC = 

Y12M_OBJ          = $(Y12M_F:.f=.$(OBJ_EXT))

#$(LIBY12M): $(Y12M_OBJ)
LIBY12M: $(Y12M_OBJ)
	mkdir -p $(TRILINOS_HOME)/external/y12m/lib/$(TRILINOS_TARGET)
	$(AR) $(ARFLAGS) $(LIB_WIN_FLAGS)$(LIBY12M) $(Y12M_OBJ)
# following ranlib commented out by rbl
#	$(RANLIB) $(LIBY12M)

#
# dependencies for 'f' files (none at this time)
#
#include ../../build/depends.y12m

clean:
	echo "junk" > dummy.$(OBJ_EXT)
	rm -f *.$(OBJ_EXT)

clobber:
	echo "junk" > dummy.$(OBJ_EXT)
	rm -f *.$(OBJ_EXT)  $(LIBY12M)

