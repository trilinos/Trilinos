#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = AZ_sync_test

AZ_sync_test_SOURCES = $(srcdir)/cc_main.cc $(srcdir)/az_sync.c

AZ_sync_test_DEPENDENCIES = $(top_builddir)/src/libaztecoo.a

AZ_sync_test_CXXFLAGS = -I$(top_srcdir)/src

AZ_sync_test_LDADD = -L$(top_builddir)/src $(FLIBS)


