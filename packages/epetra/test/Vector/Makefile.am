#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = Vector_test

Vector_test_SOURCES = $(srcdir)/cxx_main.cpp $(srcdir)/BuildTestProblems.cpp \
                  $(srcdir)/ExecuteTestProblems.cpp \
		  $(srcdir)/ExecuteTestProblems.h $(srcdir)/BuildTestProblems.h

#Vector_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

Vector_test_CXXFLAGS = -I$(top_srcdir)/src

if BUILD_BLAS
REFERENCE_BLAS_LIBS=-L$(top_builddir)/../../external/blas -lblas
endif

Vector_test_LDADD = -L$(top_builddir)/src -lepetra $(BLAS_LIBS) \
	$(REFERENCE_BLAS_LIBS) $(LIBS) $(FLIBS)
