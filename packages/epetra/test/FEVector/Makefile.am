#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = FEVector_test

FEVector_test_SOURCES = $(srcdir)/cxx_main.cpp \
                  $(srcdir)/ExecuteTestProblems.cpp \
		  $(srcdir)/ExecuteTestProblems.h 

FEVector_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

FEVector_test_CXXFLAGS = -I$(top_srcdir)/src

FEVector_test_LDADD = -L$(top_builddir)/src -lepetra $(BLAS_LIBS) \
	$(LIBS) $(FLIBS)
