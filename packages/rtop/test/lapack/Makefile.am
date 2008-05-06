
EXEEXT = .exe

include $(top_builddir)/Makefile.export.rtop

bin_PROGRAMS = testLapackWrappers

testLapackWrappers_SOURCES = \
	$(srcdir)/testLapackWrappers.cpp

testLapackWrappers_DEPENDENCIES = $(top_builddir)/src/librtop.a

testLapackWrappers_LDADD = $(RTOP_LIBS)

AM_CPPFLAGS = $(RTOP_INCLUDES)

# This is already added part of RTOP_LIBS and therefore automake does not need to add this again!
LIBS =
