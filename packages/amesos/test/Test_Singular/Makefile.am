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
TEST = Test_Singular

if HAVE_AMESOS_GALERI
if HAVE_AMESOS_KLU
noinst_PROGRAMS = $(TEST) Dummy
else
noinst_PROGRAMS = 
endif
else
noinst_PROGRAMS = 
endif

Dummy$(EXEEXT): $(srcdir)/Dummy
	if [ $(srcdir) != '.' ]; then cp  $(srcdir)/bcsstm05_ns.rua bcsstm05_ns.rua; fi

Dummy_SOURCES = Dummy

Test_Singular_SOURCES      = $(srcdir)/cxx_main.cpp
Test_Singular_DEPENDENCIES = $(top_builddir)/src/libamesos.a
Test_Singular_CXXFLAGS     = $(EXAMPLE_INCLUDES)
Test_Singular_LDADD = $(EXAMPLE_LIBS)

EXTRA_DIST = bcsstm05_ns.rua

