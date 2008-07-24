
EXEEXT = .exe

include $(top_builddir)/Makefile.export.rtop

bin_PROGRAMS = RTOpPack_Types_test baseUnitTests

RTOpPack_Types_test_SOURCES = \
	$(srcdir)/RTOpPack_Types_test.cpp
RTOpPack_Types_test_DEPENDENCIES = $(top_builddir)/src/librtop.a
RTOpPack_Types_test_LDADD = $(RTOP_LIBS)

baseUnitTests_SOURCES = \
  $(srcdir)/PrimitiveTypeTraits_UnitTests.cpp \
  $(top_srcdir)/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
baseUnitTests_DEPENDENCIES = $(top_builddir)/src/librtop.a $(top_builddir)/../teuchos/src/libteuchos.a
baseUnitTests_LDADD = $(RTOP_LIBS)

AM_CPPFLAGS = $(RTOP_INCLUDES)

# This is already added part of RTOP_LIBS and therefore automake does not need to add this again!
LIBS =
