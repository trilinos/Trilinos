#@HEADER
#@HEADER

include $(top_builddir)/Makefile.export.amesos

if USING_GNUMAKE
EXAMPLE_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(AMESOS_LIBS) $(AMESOS_EXTRA_LIBS))
EXAMPLE_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(AMESOS_INCLUDES) $(AMESOS_EXTRA_INCLUDES))
else
EXAMPLE_LIBS = $(AMESOS_LIBS) $(AMESOS_EXTRA_LIBS)
EXAMPLE_INCLUDES = $(AMESOS_INCLUDES) $(AMESOS_EXTRA_INCLUDES)
endif

EXEEXT = .exe
TEST = Test_LAPACK

if HAVE_AMESOS_LAPACK
noinst_PROGRAMS = $(TEST)
else
noinst_PROGRAMS =
endif

Test_LAPACK_SOURCES      = $(srcdir)/cxx_main.cpp
Test_LAPACK_DEPENDENCIES = $(top_builddir)/src/libamesos.a
Test_LAPACK_CXXFLAGS     = $(EXAMPLE_INCLUDES)
Test_LAPACK_LDADD        = $(EXAMPLE_LIBS)
