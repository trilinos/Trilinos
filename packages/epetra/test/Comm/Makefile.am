#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = Comm_test

Comm_test_SOURCES = $(srcdir)/cxx_main.cpp

Comm_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

Comm_test_CXXFLAGS = -I$(top_srcdir)/src
#Unsure why this test need to link to the BLAS on an IBM.
Comm_test_LDADD = -L$(top_builddir)/src -lepetra $(BLAS_LIBS) \
	$(LIBS) $(FLIBS)


