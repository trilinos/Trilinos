#@HEADER
#@HEADER
#np# This file is an example of a typical Makefile.am to build a test exe.  

EXEEXT = .exe

#np# List all executable names here
noinst_PROGRAMS = LAPACK_test

#np# One SOURCES primary is needed for EACH exe.  (Notice for all of the following 
#np# primaries that the pattern is "(Exe_name)_(PRIMARY_NAME) = ").  List all source 
#np# files needed to compile the exe.
LAPACK_test_SOURCES = $(srcdir)/cxx_main.cpp

#np# One DEPENDENCIES primary is needed for each exe.  List any libraries that the 
#np# exe depends on.  (This primary is not absolutely essential.)
LAPACK_test_DEPENDENCIES =

#np# One CXXFLAGS primary can be defined for each exe.  Include any -I's necessary for 
#np# building the exe.
LAPACK_test_CXXFLAGS = -I$(top_srcdir)/src

#np# One LDADD primary can be defined for each exe.  Include any LDFLAGS necessary for 
#np# building the exe.
LAPACK_test_LDADD = $(BLAS_LIBS) $(LIBS) $(FLIBS)


