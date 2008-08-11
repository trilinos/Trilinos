
EXEEXT = .exe

include $(top_builddir)/Makefile.export.rtop

bin_PROGRAMS = opsUnitTests

opsUnitTests_SOURCES = \
  $(srcdir)/ROpCountNanInf_UnitTests.cpp \
  $(srcdir)/ROpDotProd_UnitTests.cpp \
  $(srcdir)/ROpGetSubVector_UnitTests.cpp \
  $(srcdir)/ROpMaxIndexLessThanBound_UnitTests.cpp \
  $(srcdir)/ROpMaxIndex_UnitTests.cpp \
  $(srcdir)/ROpMax_UnitTests.cpp \
  $(srcdir)/ROpMinIndexGreaterThanBound_UnitTests.cpp \
  $(srcdir)/ROpMinIndex_UnitTests.cpp \
  $(srcdir)/ROpMin_UnitTests.cpp \
  $(srcdir)/ROpNorm1_UnitTests.cpp \
  $(srcdir)/ROpNorm2_UnitTests.cpp \
  $(srcdir)/ROpNormInf_UnitTests.cpp \
  $(srcdir)/ROpSum_UnitTests.cpp \
  $(srcdir)/ROpWeightedNorm2_UnitTests.cpp \
  $(srcdir)/TOpSetSubVector_UnitTests.cpp \
  $(srcdir)/TOpAbs_UnitTests.cpp \
  $(srcdir)/TOpAddScalar_UnitTests.cpp \
  $(srcdir)/TOpAssignScalar_UnitTests.cpp \
  $(srcdir)/TOpAssignVectors_UnitTests.cpp \
  $(srcdir)/TOpAXPY_UnitTests.cpp \
  $(srcdir)/TOpEleWiseConjProd_UnitTests.cpp \
  $(srcdir)/TOpEleWiseDivide_UnitTests.cpp \
  $(srcdir)/TOpEleWiseProd_UnitTests.cpp \
  $(srcdir)/TOpEleWiseProdUpdate_UnitTests.cpp \
  $(srcdir)/TOpLinearCombination_UnitTests.cpp \
  $(srcdir)/TOpRandomize_UnitTests.cpp \
  $(srcdir)/TOpReciprocal_UnitTests.cpp \
  $(srcdir)/TOpScaleVector_UnitTests.cpp \
  $(top_srcdir)/test/support/supportUnitTestsHelpers.cpp \
  $(top_srcdir)/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp

opsUnitTests_DEPENDENCIES = $(top_builddir)/src/librtop.a $(top_builddir)/../teuchos/src/libteuchos.a
opsUnitTests_LDADD = $(RTOP_LIBS)
AM_CPPFLAGS = $(RTOP_INCLUDES) -I$(top_srcdir)/test/support

# This is already added part of RTOP_LIBS and therefore automake does not need to add this again!
LIBS =
