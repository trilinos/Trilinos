#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = NOX_NewTest

NOX_NewTest_SOURCES = $(srcdir)/test.C

NOX_NewTest_DEPENDENCIES = $(top_builddir)/src/libnox.a \
                           $(top_builddir)/src-lapack/libnoxlapack.a

NOX_NewTest_CXXFLAGS = -I$(top_srcdir)/src \
                       -I$(top_srcdir)/src-lapack \
                       -I$(top_srcdir)/test/utils

#Unsure why this test need to link to the BLAS on an IBM.
NOX_NewTest_LDADD = -L$(top_builddir)/src -lnox $(BLAS_LIBS) \
                    -L$(top_builddir)/src-lapack -lnoxlapack \
	            $(LIBS) $(FLIBS)


