#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = NOX_Group

NOX_Group_SOURCES = $(srcdir)/test.C

NOX_Group_DEPENDENCIES = \
                         $(top_builddir)/src/libnox.a \
                         $(top_builddir)/src-epetra/libnoxepetra.a

NOX_Group_CXXFLAGS = \
                     -I$(top_srcdir)/src \
                     -I$(top_srcdir)/src-epetra \
                     -I$(top_srcdir)/test/utils \
                     $(TRILINOS_CXXFLAGS)

#Unsure why this test need to link to the BLAS on an IBM.
NOX_Group_LDADD = \
                  -L$(top_builddir)/src -lnox $(BLAS_LIBS) \
                  -L$(top_builddir)/src-epetra -lnoxepetra \
	          $(LIBS) $(FLIBS)


