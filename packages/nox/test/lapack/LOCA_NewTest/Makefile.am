#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = LOCA_NewTest

LOCA_NewTest_SOURCES = $(srcdir)/test.C

LOCA_NewTest_DEPENDENCIES = \
			  $(top_builddir)/src/libnox.a \
                          $(top_builddir)/src-lapack/libnoxlapack.a \
			  $(top_builddir)/src-loca/src/libloca.a \
                          $(top_builddir)/src-loca/src-lapack/liblocalapack.a

LOCA_NewTest_CXXFLAGS = \
                      -I$(top_srcdir)/src \
                      -I$(top_srcdir)/src-lapack \
                      -I$(top_srcdir)/src-loca/src \
                      -I$(top_srcdir)/src-loca/src-lapack \
                      -I$(top_srcdir)/test/utils

#Unsure why this test need to link to the BLAS on an IBM.
LOCA_NewTest_LDADD = \
		   -L$(top_builddir)/src -lnox $(BLAS_LIBS) \
                   -L$(top_builddir)/src-lapack -lnoxlapack \
                   -L$(top_builddir)/src-loca/src -lloca $(BLAS_LIBS) \
                   -L$(top_builddir)/src-loca/src-lapack -llocalapack \
	           $(LIBS) $(FLIBS)


