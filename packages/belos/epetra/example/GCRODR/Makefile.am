#@HEADER
#@HEADER

include $(top_builddir)/Makefile.export.belos

if USING_GNUMAKE
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(BELOS_TEST_LIBS))
EXPORT_INC_PATH = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(BELOS_TEST_INCLUDES))
else
EXPORT_LIBS = $(BELOS_TEST_LIBS)
EXPORT_INC_PATH = $(BELOS_TEST_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INC_PATH)

EXEEXT = .exe

if BUILD_TRIUTILS
if BUILD_IFPACK
BELOS_GCRODR_IFPACK_EX = prec_gcrodr_hb 
endif
if BUILD_EPETRAEXT
BELOS_GCRODR_NOPREC_EX = gcrodr_hb Dummy 
endif
endif

noinst_PROGRAMS = $(BELOS_GCRODR_NOPREC_EX) $(BELOS_GCRODR_IFPACK_EX)

Dummy$(EXEEXT): $(srcdir)/Dummy
	if [ $(srcdir) != '.' ]; then cp  $(srcdir)/sherman5.hb sherman5.hb; fi

Dummy_SOURCES = Dummy

gcrodr_hb_SOURCES = $(srcdir)/GCRODREpetraExFile.cpp $(srcdir)/createEpetraProblem.cpp

prec_gcrodr_hb_SOURCES = $(srcdir)/PrecGCRODREpetraExFile.cpp $(srcdir)/createEpetraProblem.cpp

gcrodr_hb_DEPENDENCIES = \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

prec_gcrodr_hb_DEPENDENCIES = \
			$(top_builddir)/../ifpack/src/libifpack.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

gcrodr_hb_CXXFLAGS = -I$(top_srcdir)/src 

prec_gcrodr_hb_CXXFLAGS = -I$(top_srcdir)/src 

gcrodr_hb_LDADD = $(EXPORT_LIBS)

prec_gcrodr_hb_LDADD = $(EXPORT_LIBS)

# Dummy for new export makefiles
LIBS = 

EXTRA_DIST = \
        $(srcdir)/createEpetraProblem.hpp
