
EXEEXT = .exe

include $(top_builddir)/Makefile.export.rtop

bin_PROGRAMS = supportUnitTests

supportUnitTests_SOURCES = \
  $(srcdir)/validate_apply_op_UnitTests.cpp \
  $(srcdir)/supportUnitTestsHelpers.cpp \
  $(top_srcdir)/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp

supportUnitTests_DEPENDENCIES = $(top_builddir)/src/librtop.a $(top_builddir)/../teuchos/src/libteuchos.a

supportUnitTests_LDADD = $(RTOP_LIBS)

AM_CPPFLAGS = $(RTOP_INCLUDES)

# This is already added part of RTOP_LIBS and therefore automake does not need to add this again!
LIBS =


EXTRA_DIST = $(srcdir)/supportUnitTestsHelpers.hpp
