#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = LOCA_NewTest

LOCA_NewTest_SOURCES = $(srcdir)/test.C

LOCA_NewTest_DEPENDENCIES = \
                          $(top_builddir)/src/libnox.a \
                          $(top_builddir)/src-epetra/libnoxepetra.a \
                          $(top_builddir)/src-loca/src/libloca.a \
                          $(top_builddir)/src-loca/src-epetra/liblocaepetra.a

LOCA_NewTest_CXXFLAGS = \
                      -I$(top_srcdir)/src \
                      -I$(top_srcdir)/src-epetra \
                      -I$(top_srcdir)/src-loca/src \
                      -I$(top_srcdir)/src-loca/src-epetra \
                      -I$(top_srcdir)/test/utils \
                      $(TRILINOS_CXXFLAGS)

#Unsure why this test need to link to the BLAS on an IBM.
LOCA_NewTest_LDADD = \
                   -L$(top_builddir)/src -lnox $(BLAS_LIBS) \
                   -L$(top_builddir)/src-epetra -lnoxepetra \
                   -L$(top_builddir)/src-loca/src -lloca $(BLAS_LIBS) \
                   -L$(top_builddir)/src-loca/src-epetra -llocaepetra \
	           $(LIBS) $(FLIBS)


