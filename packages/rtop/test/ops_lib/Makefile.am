
EXEEXT = .exe

include $(top_builddir)/Makefile.export.rtop

bin_PROGRAMS = opsUnitTests

opsUnitTests_SOURCES = \
  $(srcdir)/ROpSum_UnitTests.cpp \
  $(srcdir)/ROpDotProd_UnitTests.cpp \
  $(srcdir)/ROpNorm1_UnitTests.cpp \
  $(srcdir)/ROpNorm2_UnitTests.cpp \
  $(srcdir)/ROpNormInf_UnitTests.cpp \
  $(srcdir)/opsUnitTestsHelpers.cpp \
  $(top_srcdir)/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp

opsUnitTests_DEPENDENCIES = $(top_builddir)/src/librtop.a $(top_builddir)/../teuchos/src/libteuchos.a

opsUnitTests_LDADD = $(RTOP_LIBS)

AM_CPPFLAGS = $(RTOP_INCLUDES)

# This is already added part of RTOP_LIBS and therefore automake does not need to add this again!
LIBS =

EXTRA_DIST = opsUnitTestsHelpers.hpp
