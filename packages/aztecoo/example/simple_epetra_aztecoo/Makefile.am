#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = epetra_crsmatrix

epetra_crsmatrix_SOURCES = $(srcdir)/cxx_main.cpp

#epetra_crsmatrix_DEPENDENCIES = $(top_builddir)/src/libaztecoo.a $(top_builddir)/../epetra/src/libepetra.a

epetra_crsmatrix_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/../epetra/src \
	-I$(top_builddir)/../epetra/src

if BUILD_Y12M
REFERENCE_Y12M_LIBS=-L$(top_builddir)/../../external/y12m -ly12m
endif

if BUILD_OLDIFPACK
OLDIFPACK_LIBS=-L$(top_builddir)/../ifpack/src -lifpack
endif 

epetra_crsmatrix_LDADD = -L$(top_builddir)/src -laztecoo $(OLDIFPACK_LIBS) \
	-L$(top_builddir)/../epetra/src -lepetra $(REFERENCE_Y12M_LIBS) \
	$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)


