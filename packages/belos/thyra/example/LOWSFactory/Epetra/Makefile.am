include $(top_builddir)/Makefile.export.belos

if USING_GNUMAKE
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(BELOS_TEST_LIBS))
EXPORT_INC_PATH = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(BELOS_TEST_INCLUDES))
else
EXPORT_LIBS = $(BELOS_TEST_LIBS)
EXPORT_INC_PATH = $(BELOS_TEST_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INC_PATH)

EXEEXT = .exe

#
# Actual programs
#

if BUILD_EPETRAEXT
BELOS_THYRA_SOLVER_DRIVER = belos_epetra_thyra_lowsf_hb 
endif

noinst_PROGRAMS = $(BELOS_THYRA_SOLVER_DRIVER) 

belos_epetra_thyra_lowsf_hb_SOURCES = $(srcdir)/belos_epetra_thyra_lowsf_hb.cpp 
belos_epetra_thyra_lowsf_hb_DEPENDENCIES= $(common_dependencies)
belos_epetra_thyra_lowsf_hb_LDADD = $(common_ldadd)

Example.mtx:
	cp  $(top_srcdir)/../belos/test/BlockGmres/orsirr1.hb orsirr1.hb                               
	cp  $(top_srcdir)/../belos/test/BlockCG/bcsstk14.hb bcsstk14.hb                              

#
# Common stuff
#

common_dependencies = $(top_builddir)/src/libbelos.a Example.mtx

common_ldadd = $(EXPORT_LIBS)

# This is already added part of BELOS_LIBS and therefore automake does not need to add this again!
LIBS =
