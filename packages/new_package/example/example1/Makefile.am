#@HEADER
#@HEADER
#np# This file is an example of a typical Makefile.am to build an example exe. 

EXEEXT = .exe

#np# List the names of all exe's in this directory here
noinst_PROGRAMS = example

#np# One SOURCES primary is needed for EACH exe listed above.  (Notice for all 
#np# of the  following primaries that the pattern is 
#np# "(Exe_name)_(PRIMARY_NAME) = ").
#np# List all source files needed to compile the exe.
example_SOURCES = $(srcdir)/cxx_main.cpp

#np# One DEPENDENCIES primary is needed for each exe listed above.  List any
#np# libraries that the exe depends on.  (This primary is not absolutely
#np# essential.)
example_DEPENDENCIES = $(top_builddir)/src/libepetra.a

#np# One CXXFLAGS primary can be defined for each exe listed above.  Include 
#np# any -I's necessary for building the exe.
example_CXXFLAGS = -I$(top_srcdir)/src

#np# One LDADD primary can be defined for each exe listed above.  Include any
#np# LDFLAGS  necessary for building the exe.
example_LDADD = $(top_builddir)/src/libepetra.a $(BLAS_LIBS) \
	$(LIBS) $(FLIBS)
