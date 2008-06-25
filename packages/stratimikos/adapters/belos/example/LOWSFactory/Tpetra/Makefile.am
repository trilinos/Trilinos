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

if ENABLE_STRATIMIKOS_TRIUTILS
BELOS_THYRA_SOLVER_DRIVER_HB = belos_tpetra_thyra_lowsf_hb 
endif

# noinst_PROGRAMS = belos_tpetra_thyra_lowsf_lap \
#                   belos_tpetra_thyra_lowsf_prec_lap \
#                   $(BELOS_THYRA_SOLVER_DRIVER_HB) 
noinst_PROGRAMS = 
                 
                

belos_tpetra_thyra_lowsf_lap_SOURCES = $(srcdir)/belos_tpetra_thyra_lowsf_lap.cpp 
belos_tpetra_thyra_lowsf_lap_DEPENDENCIES= $(common_dependencies)
belos_tpetra_thyra_lowsf_lap_LDADD = $(common_ldadd)

belos_tpetra_thyra_lowsf_prec_lap_SOURCES = $(srcdir)/belos_tpetra_thyra_lowsf_prec_lap.cpp 
belos_tpetra_thyra_lowsf_prec_lap_DEPENDENCIES= $(common_dependencies)
belos_tpetra_thyra_lowsf_prec_lap_LDADD = $(common_ldadd)

belos_tpetra_thyra_lowsf_hb_SOURCES = $(srcdir)/belos_tpetra_thyra_lowsf_hb.cpp 
belos_tpetra_thyra_lowsf_hb_DEPENDENCIES= $(common_dependencies)
belos_tpetra_thyra_lowsf_hb_LDADD = $(common_ldadd)

mhd1280b.cua:
	cp  $(top_srcdir)/../belos/test/MVOPTester/mhd1280b.cua mhd1280b.cua

#
# Common stuff
#

common_dependencies = $(top_builddir)/../belos/src/libbelos.a mhd1280b.cua

common_ldadd = $(EXPORT_LIBS)

# This is already added part of BELOS_LIBS and therefore automake does not
# need to add this again!
LIBS = 

EXTRA_DIST = \
	$(srcdir)/MyOperator.hpp
