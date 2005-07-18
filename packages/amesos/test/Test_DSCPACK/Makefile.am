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

AM_CPPFLAGS = $(EXAMPLE_INCLUDES)

EXEEXT = .exe
TEST = Test_DSCPACK

noinst_PROGRAMS = $(TEST)

Test_DSCPACK_SOURCES = $(srcdir)/cxx_main.cpp
Test_DSCPACK_DEPENDENCIES =  all-libs
Test_DSCPACK_LDADD = $(EXAMPLE_LIBS)

all-libs:
	cd $(top_builddir)/src ; $(MAKE) all-libs 
