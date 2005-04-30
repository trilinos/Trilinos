##
## Example external makefile that shows how to compile and
## link against an installed version of Thyra (part of Trilinos)
##

# Points to where Trilinos is installed
TRILINOS_INSTALL_DIR = ../..

# Include the makefile fragment that specify the compile/link commands 
include $(TRILINOS_INSTALL_DIR)/include/Makefile.export.thyra.macros
# Include the makefile fragement that defines the include paths and libraries
include $(TRILINOS_INSTALL_DIR)/include/Makefile.export.thyra

#
# Get the compile and link flags
#

# For platforms where perl and GNUMake are used
COMPILE_FLAGS = $(shell perl $(TRILINOS_INSTALL_DIR)/include/strip_dup_incl_paths.pl \
	$(THYRA_CXXFLAGS) $(THYRA_DEFS) $(THYRA_CPPFLAGS) $(THYRA_INCLUDES))
# For platforms where perl and/or GNUMake are not available
#COMPILE_FLAGS = $(THYRA_CXXFLAGS) $(THYRA_DEFS) $(THYRA_CPPFLAGS) $(THYRA_INCLUDES)
# For platforms where perl and GNUMake are used
LINK_FLAGS = $(shell perl $(TRILINOS_INSTALL_DIR)/include/strip_dup_libs.pl $(THYRA_LIBS))
# For platforms where perl and/or GNUMake are not available
#LINK_FLAGS = $(THYRA_LIBS)

sillyCgSolve_serial.exe : sillyCgSolve_serial.o
	$(THYRA_CXXLD) $(THYRA_CXXFLAGS) sillyCgSolve_serial.o -o sillyCgSolve_serial.exe $(LINK_FLAGS)

sillyCgSolve_serial.o : sillyCgSolve_serial.cpp
	$(THYRA_CXX) $(COMPILE_FLAGS) -c sillyCgSolve_serial.cpp

clean:
	rm -f *.o *.exe *~
