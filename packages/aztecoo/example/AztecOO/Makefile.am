#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = AztecOO_test

AztecOO_test_SOURCES = $(srcdir)/cc_main.cc $(srcdir)/az_sync.c

AztecOO_test_DEPENDENCIES = $(top_builddir)/src/libaztecoo.a

AztecOO_test_CXXFLAGS = -I$(top_srcdir)/src

AztecOO_test_LDADD = -L$(top_builddir)/src $(FLIBS)


