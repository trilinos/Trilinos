#@HEADER
#@HEADER

EXEEXT = .exe

noinst_PROGRAMS = amd_simple

amd_simple_CFLAGS = -I$(top_srcdir)/amd/src

amd_simple_SOURCES = $(srcdir)/amd_simple.c

amd_simple_DEPENDENCIES = \
	$(top_builddir)/amd/src/libamd.a 

amd_simple_LDADD = $(top_builddir)/amd/src/libamd.a $(FLIBS)


