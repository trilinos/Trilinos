#@HEADER
#@HEADER

EXEEXT = .exe

INTERFACE = $(srcdir)/Interface.C

noinst_PROGRAMS = testNonlinear

testNonlinear_SOURCES = \
                        $(INTERFACE) \
                        $(srcdir)/testNonlinear.C

testNonlinear_DEPENDENCIES = \
                             $(top_builddir)/src/libnox.a \
                             $(top_builddir)/src-epetra/libnoxepetra.a

testNonlinear_CXXFLAGS = \
                         -I$(top_srcdir)/src \
                         -I$(top_srcdir)/src-epetra \
                         -I$(top_srcdir)/test/utils \
                         $(TRILINOS_CXXFLAGS)

testNonlinear_LDADD = \
                      -L$(top_builddir)/src -lnox $(BLAS_LIBS) \
                      -L$(top_builddir)/src-epetra -lnoxepetra \
	              $(LIBS) $(FLIBS)


