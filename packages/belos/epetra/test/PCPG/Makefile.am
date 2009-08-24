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
if BUILD_EPETRAEXT
if BUILD_ML
BELOS_BCG_ML_EX = pcpg_epetraex
endif
endif
endif

noinst_PROGRAMS = $(BELOS_BCG_ML_EX)

pcpg_epetraex_SOURCES = $(srcdir)/test_pcpg_epetraex.cpp

pcpg_epetraex_DEPENDENCIES = \
			$(top_builddir)/../ml/src/libml.a \
			$(top_builddir)/../epetraext/src/libepetraext.a \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

pcpg_epetraexprec_CXXFLAGS = -I$(top_srcdir)/src 

pcpg_epetraex_LDADD = $(EXPORT_LIBS)

# Dummy for new export makefiles
LIBS = 

