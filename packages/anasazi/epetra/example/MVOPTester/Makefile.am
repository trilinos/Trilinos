#@HEADER
#@HEADER

#np# Specify locations of include files within Trilinos, but outside np
#np# Add packages that np is dependent on; remove epetra if np does not need it

include $(top_builddir)/Makefile.export.anasazi

if USING_GNUMAKE
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(ANASAZI_TEST_LIBS))
EXPORT_INC_PATH = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(ANASAZI_TEST_INCLUDES))
else
EXPORT_LIBS = $(ANASAZI_TEST_LIBS)
EXPORT_INC_PATH = $(ANASAZI_TEST_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INC_PATH)

EXEEXT = .exe

noinst_PROGRAMS = MVOPTesterEx

MVOPTesterEx_SOURCES = $(srcdir)/MVOPTesterEx.cpp

#np# Specify locations of include files within Trilinos, but outside np
#np# Add packages that np is dependent on; remove epetra if np does not need it

#np# - List libraries that np is dependent on.  The DEPENDENCIES primary is
#np# optional.
MVOPTesterEx_DEPENDENCIES = \
        $(top_builddir)/src/libanasazi.a \
        $(top_builddir)/../epetra/src/libepetra.a \
        $(top_builddir)/../teuchos/src/libteuchos.a 

MVOPTesterEx_CXXFLAGS = \
			-I$(top_srcdir)/src \
			-I$(top_srcdir)/epetra/util/ModeLaplace

#np# List libraries that np is dependent on.
MVOPTesterEx_LDADD = \
        $(EXPORT_LIBS)

#np# Set "SUBDIRS= (list all subdirectories that are autotool'ed here)"
#np# List all subdirectories even if one or more are configured/built only
#np# conditionally.
#np# For packages whose examples are located in this directory, rather than
#np# in a subdirectory of this directory, use the Makefile.am located in
#np# example1/ rather than this file.
# SUBDIRS = example1
SUBDIRS = 

# Dummy for new export makefiles
LIBS = 
