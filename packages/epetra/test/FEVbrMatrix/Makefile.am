#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = FEVbrMatrix

FEVbrMatrix_test_SOURCES = $(srcdir)/cxx_main.cpp \
                  $(srcdir)/ExecuteTestProblems.cpp \
		  $(srcdir)/ExecuteTestProblems.h 

FEVbrMatrix_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

FEVbrMatrix_test_CXXFLAGS = -I$(top_srcdir)/src

FEVbrMatrix_test_LDADD = -L$(top_builddir)/src -lepetra $(BLAS_LIBS) \
	$(LIBS) $(FLIBS)
