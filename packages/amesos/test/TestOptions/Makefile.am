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

TEST = TestOptions 

noinst_PROGRAMS = $(TEST) 

# if HAVE_AMESOS_EPETRAEXT
#    EPETRA_EXT_LIB = 	$(top_builddir)/../epetraext/src/libepetraext.a 
#    EPETRA_EXT_INCLUDE =	-I$(top_srcdir)/../epetraext/src/transform -I$(top_builddir)/../epetraext/src/transform -I$(top_srcdir)/../epetraext/src -I$(top_builddir)/../epetraext/src 
# else
#    EPETRA_EXT_LIB = 
#    EPETRA_EXT_INCLUDE = 
# endif 

EXTRA_AMESOS_TRIUTILS_TEST = $(srcdir)/TestOptions.cpp  \
        $(srcdir)/CreateTridi.cpp \
        $(srcdir)/CreateTridi.h \
        $(srcdir)/TestSuperludist.cpp \
        $(srcdir)/TestScalapack.cpp \
        $(srcdir)/TestScalapack.h \
        $(srcdir)/TestKlu.cpp \
        $(srcdir)/TestKlu.h \
        $(srcdir)/TestSuperludist.h \
        $(srcdir)/NewMatNewMap.cpp \
        $(srcdir)/NewMatNewMap.h \
        $(srcdir)/TestAllClasses.cpp \
        $(srcdir)/TestAllClasses.h \
        $(srcdir)/PartialFactorization.cpp \
        $(srcdir)/PartialFactorization.h \
        $(srcdir)/PerformOneSolveAndTest.cpp \
        $(srcdir)/PerformOneSolveAndTest.h \
        $(srcdir)/TestOtherClasses.cpp \
        $(srcdir)/TestOtherClasses.h

EXTRA_TRIVIAL_TEST = $(srcdir)/DoNothing.cpp

if HAVE_AMESOS_TRIUTILS
TestOptions_SOURCES = $(EXTRA_AMESOS_TRIUTILS_TEST)
else
TestOptions_SOURCES = $(EXTRA_TRIVIAL_TEST)
endif

EXTRA_TestOptions_SOURCES = \
	$(EXTRA_AMESOS_TRIUTILS_TEST) \
	$(EXTRA_TRIVIAL_TEST)

TestOptions_DEPENDENCIES = $(top_builddir)/src/libamesos.a
TestOptions_LDADD = $(EXAMPLE_LIBS)
