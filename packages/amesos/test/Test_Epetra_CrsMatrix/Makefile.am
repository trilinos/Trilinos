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
TEST = Test_Epetra_CrsMatrix

noinst_PROGRAMS = $(TEST)

Test_Epetra_CrsMatrix_SOURCES      = $(srcdir)/cxx_main.cpp
Test_Epetra_CrsMatrix_DEPENDENCIES = $(top_builddir)/src/libamesos.a
Test_Epetra_CrsMatrix_CXXFLAGS     = $(EXAMPLE_INCLUDES)
Test_Epetra_CrsMatrix_LDADD        = $(EXAMPLE_LIBS)
