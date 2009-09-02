
include $(top_builddir)/Makefile.export.stratimikos

_EXPORT_INC_PATH = $(STRATIMIKOS_TEST_INCLUDES)
_EXPORT_LIBS = $(STRATIMIKOS_TEST_LIBS)

if USING_GNUMAKE
EXPORT_INC_PATH = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(_EXPORT_INC_PATH))
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(_EXPORT_LIBS))
else
EXPORT_INC_PATH = $(_EXPORT_INC_PATH)
EXPORT_LIBS = $(_EXPORT_LIBS)
endif

AM_CPPFLAGS = $(EXPORT_INC_PATH)

EXEEXT = .exe

#np# List all executable names here
noinst_PROGRAMS = MVOPTesterThyra_test

#np# One SOURCES primary is needed for EACH exe.  (Notice for all of the following 
#np# primaries that the pattern is "(Exe_name)_(PRIMARY_NAME) = ").  List all source 
#np# files needed to compile the exe.
MVOPTesterThyra_test_SOURCES = $(srcdir)/cxx_main.cpp

#np# One DEPENDENCIES primary is needed for each exe.  List any libraries that the 
#np# exe depends on.  (This primary is not absolutely essential.)
MVOPTesterThyra_test_DEPENDENCIES = \
	$(top_builddir)/../belos/src/libbelos.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(top_builddir)/../teuchos/src/libteuchos.a

#np# One CXXFLAGS primary can be defined for each exe.  Include any -I's necessary for 
#np# building the exe.
MVOPTesterThyra_test_CXXFLAGS = \
	-I$(top_srcdir)/thyra/src \
	-I$(top_srcdir)/src

#np# One LDADD primary can be defined for each exe.  Include any LDFLAGS necessary for 
#np# building the exe.
MVOPTesterThyra_test_LDADD = $(EXPORT_LIBS) 

# This is already added part of BELOS_LIBS and therefore automake does not
# need to add this again!
LIBS = 


