#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = FECrsMatrix_test

FECrsMatrix_test_SOURCES = $(srcdir)/cxx_main.cpp \
                  $(srcdir)/ExecuteTestProblems.cpp \
		  $(srcdir)/ExecuteTestProblems.h 

FECrsMatrix_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

FECrsMatrix_test_CXXFLAGS = -I$(top_srcdir)/src

FECrsMatrix_test_LDADD = -L$(top_builddir)/src -lepetra $(BLAS_LIBS) \
	$(LIBS) $(FLIBS)
