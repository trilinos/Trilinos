include $(top_builddir)/Makefile.export.stratimikos

_EXPORT_INC_PATH = $(STRATIMIKOS_TEST_INCLUDES)
_EXPORT_LIBS = $(STRATIMIKOS_TEST_LIBS)

if USING_GNUMAKE
EXPORT_INC_PATH = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(_EXPORT_INC_PATH))
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(_EXPORT_LIBS))
else
EXPORT_INC_PATH = $(_EXPORT_INC_PATH)
EXPORT_LIBS = $(_EXPORT_LIBS)
endif

AM_CPPFLAGS = $(EXPORT_INC_PATH)

EXEEXT = .exe

#
# Actual programs
#

if ENABLE_STRATIMIKOS_THYRA_EPETRAEXT
BELOS_THYRA_SOLVER_DRIVER = belos_epetra_thyra_lowsf_hb 
endif

noinst_PROGRAMS = $(BELOS_THYRA_SOLVER_DRIVER) 

belos_epetra_thyra_lowsf_hb_SOURCES = $(srcdir)/belos_epetra_thyra_lowsf_hb.cpp 
belos_epetra_thyra_lowsf_hb_DEPENDENCIES= $(common_dependencies)
belos_epetra_thyra_lowsf_hb_LDADD = $(common_ldadd)

orsirr1.hb:
	cp  $(top_srcdir)/../belos/example/BlockGmres/orsirr1.hb orsirr1.hb                               
	cp  $(top_srcdir)/../belos/example/BlockCG/bcsstk14.hb bcsstk14.hb

#
# Common stuff
#

common_dependencies = $(top_builddir)/../belos/src/libbelos.a orsirr1.hb

common_ldadd = $(EXPORT_LIBS)

# This is already added part of BELOS_LIBS and therefore automake does not
# need to add this again!
LIBS =
