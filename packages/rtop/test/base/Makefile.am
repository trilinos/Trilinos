
EXEEXT = .exe

include $(top_builddir)/Makefile.export.rtop

bin_PROGRAMS = RTOpPack_Types_test

RTOpPack_Types_test_SOURCES = \
	$(srcdir)/RTOpPack_Types_test.cpp

RTOpPack_Types_test_DEPENDENCIES = $(top_builddir)/src/librtop.a

RTOpPack_Types_test_LDADD = $(RTOP_LIBS)

AM_CPPFLAGS = $(RTOP_INCLUDES)

# This is already added part of RTOP_LIBS and therefore automake does not need to add this again!
LIBS =
