
EXEEXT = .exe

include $(top_builddir)/Makefile.export.rtop

bin_PROGRAMS = runSpmdTests

runSpmdTests_SOURCES = \
	$(srcdir)/runSpmdTests.cpp

runSpmdTests_DEPENDENCIES = $(top_builddir)/src/librtop.a

runSpmdTests_LDADD = $(RTOP_LIBS)

AM_CPPFLAGS = $(RTOP_INCLUDES)

# This is already added part of RTOP_LIBS and therefore automake does not need to add this again!
LIBS =
