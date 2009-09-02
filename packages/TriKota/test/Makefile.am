
EXEEXT = .exe

all-local: dakota.in

dakota.in:
	cp  $(srcdir)/_dakota.in dakota.in

#
# Actual programs
#

noinst_PROGRAMS = Simple

Simple_SOURCES = \
  Simple_Main.cpp \
  Simple_ModelEval.cpp
Simple_DEPENDENCIES = $(common_dependencies)
Simple_LDADD = $(common_ldadd)

#
# Common stuff
#

common_dependencies = \
  $(top_builddir)/src/libtrikota.a

include $(top_builddir)/Makefile.export.TriKota

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(TRIKOTA_INCLUDES))
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(TRIKOTA_LIBS))
else
EXPORT_INCLUDES = $(TRIKOTA_INCLUDES)
EXPORT_LIBS = $(TRIKOTA_LIBS)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

common_ldadd = $(EXPORT_LIBS)

# This is already added part of TRIKOTA_LIBS and therefore automake does not need to add this again!
#LIBS =

EXTRA_DIST = definition Simple_ModelEval.hpp _dakota.in
