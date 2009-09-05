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
BELOS_RCG_NOPREC_EX = rcg_hb Dummy 
endif

noinst_PROGRAMS = $(BELOS_RCG_NOPREC_EX)

Dummy$(EXEEXT): $(srcdir)/Dummy
	cp $(top_srcdir)/epetra/example/RCG/bcsstk14.hb bcsstk14.hb

Dummy_SOURCES = Dummy

rcg_hb_SOURCES = $(srcdir)/test_rcg_hb.cpp $(srcdir)/createEpetraProblem.cpp

rcg_hb_DEPENDENCIES = \
			$(top_builddir)/../epetra/src/libepetra.a \
			$(top_builddir)/../teuchos/src/libteuchos.a

rcg_hb_CXXFLAGS = -I$(top_srcdir)/src 

rcg_hb_LDADD = $(EXPORT_LIBS)

# Dummy for new export makefiles
LIBS = 

EXTRA_DIST = \
        $(srcdir)/createEpetraProblem.hpp
