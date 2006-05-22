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

if BUILD_TRIUTILS
BELOS_THYRA_SOLVER_DRIVER_HB = belos_tpetra_thyra_lowsf_hb 
endif

noinst_PROGRAMS = belos_tpetra_thyra_lowsf_lap \
                  belos_tpetra_thyra_lowsf_prec_lap \
                  $(BELOS_THYRA_SOLVER_DRIVER_HB) 

belos_tpetra_thyra_lowsf_lap_SOURCES = $(srcdir)/belos_tpetra_thyra_lowsf_lap.cpp 
belos_tpetra_thyra_lowsf_lap_DEPENDENCIES= $(common_dependencies)
belos_tpetra_thyra_lowsf_lap_LDADD = $(common_ldadd)

belos_tpetra_thyra_lowsf_prec_lap_SOURCES = $(srcdir)/belos_tpetra_thyra_lowsf_prec_lap.cpp 
belos_tpetra_thyra_lowsf_prec_lap_DEPENDENCIES= $(common_dependencies)
belos_tpetra_thyra_lowsf_prec_lap_LDADD = $(common_ldadd)

belos_tpetra_thyra_lowsf_hb_SOURCES = $(srcdir)/belos_tpetra_thyra_lowsf_hb.cpp 
belos_tpetra_thyra_lowsf_hb_DEPENDENCIES= $(common_dependencies)
belos_tpetra_thyra_lowsf_hb_LDADD = $(common_ldadd)

Example.mtx:
	cp  $(top_srcdir)/../belos/test/BlockCG/mhd1280b.cua mhd1280b.cua

#
# Common stuff
#

common_dependencies = $(top_builddir)/src/libbelos.a Example.mtx

common_ldadd = $(EXPORT_LIBS)

# This is already added part of BELOS_LIBS and therefore automake does not need to add this again!
LIBS =
