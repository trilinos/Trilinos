
EXEEXT = .exe

#
# Actual programs
#

noinst_PROGRAMS = simpleEpetraExteranCInterface

simpleEpetraExteranCInterface_SOURCES = \
  $(srcdir)/simpleEpetraExteranCInterface.c \
  $(srcdir)/FEpetra_Map.cpp \
  $(srcdir)/FEpetra_Vector.cpp
simpleEpetraExteranCInterface_DEPENDENCIES = $(common_dependencies)
simpleEpetraExteranCInterface_LDADD = $(common_ldadd)


#
# Common stuff
#

common_dependencies = \
  $(top_builddir)/../epetra/src/libepetra.a

include $(top_builddir)/Makefile.export.epetraext

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(EPETRAEXT_INCLUDES))
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(EPETRAEXT_LIBS))
else
EXPORT_INCLUDES = $(EPETRAEXT_INCLUDES)
EXPORT_LIBS = $(EPETRAEXT_LIBS)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

common_ldadd = $(EXPORT_LIBS)

# This is already added part of EPETRAEXT_LIBS and therefore automake does not need to add this again!
LIBS =
