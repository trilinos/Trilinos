#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = NOX_Group

NOX_Group_SOURCES = $(srcdir)/test.C

NOX_Group_DEPENDENCIES = $(top_builddir)/src/libnox.a \
                         $(top_builddir)/src-lapack/libnoxlapack.a

NOX_Group_CXXFLAGS = -I$(top_srcdir)/src \
                     -I$(top_srcdir)/src-lapack \
                     -I$(top_srcdir)/test/utils

#Unsure why this test need to link to the BLAS on an IBM.
NOX_Group_LDADD = -L$(top_builddir)/src -lnox $(BLAS_LIBS) \
                  -L$(top_builddir)/src-lapack -lnoxlapack \
	          $(LIBS) $(FLIBS)


