#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = NOX_NewTest

NOX_NewTest_SOURCES = $(srcdir)/test.C

NOX_NewTest_DEPENDENCIES = \
                         $(top_builddir)/src/libnox.a \
                         $(top_builddir)/src-epetra/libnoxepetra.a

NOX_NewTest_CXXFLAGS = \
                     -I$(top_srcdir)/src \
                     -I$(top_srcdir)/src-epetra \
                     -I$(top_srcdir)/test/utils \
                     $(TRILINOS_CXXFLAGS)

#Unsure why this test need to link to the BLAS on an IBM.
NOX_NewTest_LDADD = \
                  -L$(top_builddir)/src -lnox $(BLAS_LIBS) \
                  -L$(top_builddir)/src-epetra -lnoxepetra \
	          $(LIBS) $(FLIBS)


