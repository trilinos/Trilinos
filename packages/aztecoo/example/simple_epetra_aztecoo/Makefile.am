#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = epetra_crsmatrix

epetra_crsmatrix_SOURCES = $(srcdir)/cxx_main.cpp

#epetra_crsmatrix_DEPENDENCIES = $(top_builddir)/src/libaztecoo.a $(top_builddir)/../epetra/src/libepetra.a

epetra_crsmatrix_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/../epetra/src \
	-I$(top_builddir)/../epetra/src

epetra_crsmatrix_LDADD = $(top_builddir)/src/libaztecoo.a \
	$(top_builddir)/../epetra/src/libepetra.a \
	$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)


