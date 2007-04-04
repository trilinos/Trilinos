
EXEEXT = .exe

#
# Actual programs
#

SILLY_POWER_METHOD_TPETRA = sillyPowerMethod_tpetra

noinst_PROGRAMS = $(SILLY_POWER_METHOD_TPETRA)

sillyPowerMethod_tpetra_SOURCES = $(srcdir)/sillyPowerMethod_tpetra.cpp $(common_sources)
sillyPowerMethod_tpetra_DEPENDENCIES = $(common_dependencies)
sillyPowerMethod_tpetra_LDADD = $(common_ldadd)

#
# Common stuff
#

common_sources =

common_dependencies = $(top_builddir)/src/libtpetra.a

include $(top_builddir)/Makefile.export.tpetra

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(TPETRA_INCLUDES))
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(TPETRA_LIBS))
else
EXPORT_INCLUDES = $(TPETRA_INCLUDES)
EXPORT_LIBS = $(TPETRA_LIBS)
endif

AM_CPPFLAGS = -I$(top_srcdir)/../thyra/example/operator_vector $(EXPORT_INCLUDES)

common_ldadd = $(EXPORT_LIBS)

EXTRA_DIST = $(srcdir)/createExampleTridiagTpetraLinearOp.hpp

# This is already added part of TPETRA_LIBS and therefore automake does not need to add this again!
LIBS =
