#####################################################################
## Example Makefile for a user application that does not use autoconf
## - Uses lapack concrete instations for group and vector
## - Must use gnu-make (gmake) if the "shell" command is invoked
#####################################################################

##
## Set the Trilinos install directory
##
TRILINOS_INSTALL_DIR = /home/rppawlo/trilinos_local

##
## Include any direct Trilinos library dependencies - in this case only nox
##
include $(TRILINOS_INSTALL_DIR)/include/Makefile.export.nox.macros
include $(TRILINOS_INSTALL_DIR)/include/Makefile.export.nox

##
## Use one of the following lines (2nd line is for non-gnumake platforms)
##
COMPILE_FLAGS = $(shell perl $(TRILINOS_INSTALL_DIR)/include/strip_dup_incl_paths.pl $(NOX_CXXFLAGS) $(NOX_DEFS) $(NOX_CPPFLAGS) $(NOX_INCLUDES))

#COMPILE_FLAGS = $(NOX_CXXFLAGS) $(NOX_DEFS) $(NOX_CPPFLAGS) $(NOX_INCLUDES)

##
## Use one of the following lines (2nd line is for non-gnumake platforms)
##
LINK_FLAGS = $(shell perl $(TRILINOS_INSTALL_DIR)/include/strip_dup_libs.pl $(NOX_LIBS))

#LINK_FLAGS = $(NOX_LIBS)

##
## Build your application code
##

Rosenbrock.exe: Rosenbrock.o
    $(NOX_CXXLD) $(NOX_CXXFLAGS) -o Rosenbrock.exe Rosenbrock.o $(LINK_FLAGS)

Rosenbrock.o:
    $(NOX_CXX) $(COMPILE_FLAGS) -c Rosenbrock.C

clean:
    rm -f *.o Rosenbrock.exe *~
