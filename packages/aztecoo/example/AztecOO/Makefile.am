#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = AztecOO_msr adapt_msr AztecOO_vbr

AztecOO_msr_SOURCES = $(srcdir)/cxx_main.cpp

#AztecOO_msr_DEPENDENCIES = $(top_builddir)/src/libaztecoo.a \
#	$(top_builddir)/../epetra/src/libepetra.a\
#	$(top_builddir)/../triutils/src/libtriutils.a

AztecOO_msr_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/../epetra/src \
	-I$(top_srcdir)/../triutils/src -I$(top_builddir)/../epetra/src \
	-I$(top_builddir)/../triutils/src

AztecOO_msr_LDADD = $(top_builddir)/src/libaztecoo.a\
	$(top_builddir)/../triutils/src/libtriutils.a\
	$(top_builddir)/../epetra/src/libepetra.a\
	$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

AztecOO_vbr_SOURCES = $(srcdir)/vbr_main.cpp

#AztecOO_vbr_DEPENDENCIES = $(top_builddir)/src/libaztecoo.a \
#	$(top_builddir)/../epetra/src/libepetra.a\
#	$(top_builddir)/../triutils/src/libtriutils.a

AztecOO_vbr_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/../epetra/src \
	-I$(top_srcdir)/../triutils/src -I$(top_builddir)/../epetra/src \
	-I$(top_builddir)/../triutils/src

AztecOO_vbr_LDADD = $(top_builddir)/src/libaztecoo.a\
	$(top_builddir)/../triutils/src/libtriutils.a\
	$(top_builddir)/../epetra/src/libepetra.a\
	$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)

adapt_msr_SOURCES = $(srcdir)/adapt_main.cpp

#adapt_msr_DEPENDENCIES = $(top_builddir)/src/libaztecoo.a \
#	$(top_builddir)/../epetra/src/libepetra.a\
#	$(top_builddir)/../triutils/src/libtriutils.a

adapt_msr_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/../epetra/src \
	-I$(top_srcdir)/../triutils/src -I$(top_builddir)/../epetra/src \
	-I$(top_builddir)/../triutils/src

adapt_msr_LDADD = $(top_builddir)/src/libaztecoo.a\
	$(top_builddir)/../triutils/src/libtriutils.a\
	$(top_builddir)/../epetra/src/libepetra.a\
	$(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)


