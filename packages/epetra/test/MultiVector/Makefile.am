#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = MultiVector_test

MultiVector_test_SOURCES = $(srcdir)/cxx_main.cpp \
           $(srcdir)/BuildTestProblems.cpp $(srcdir)/ExecuteTestProblems.cpp \
	   $(srcdir)/BuildTestProblems.h $(srcdir)/ExecuteTestProblems.h

MultiVector_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

MultiVector_test_CXXFLAGS = -I$(top_srcdir)/src

MultiVector_test_LDADD = -L$(top_builddir)/src -lepetra $(BLAS_LIBS) $(LIBS) $(FLIBS)
