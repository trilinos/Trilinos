# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
include $(top_builddir)/Makefile.export.ml

TEST_LIBS     = $(ML_LIBS)     $(ML_EXTRA_LIBS)
TEST_INCLUDES = $(ML_INCLUDES) $(ML_EXTRA_INCLUDES)

EXEEXT = .exe

noinst_PROGRAMS = MLAPILinearCombination

MLAPILinearCombination_SOURCES      = $(srcdir)/cxx_main.cpp
MLAPILinearCombination_DEPENDENCIES = $(top_builddir)/src/libml.a
MLAPILinearCombination_CXXFLAGS     = ${TEST_INCLUDES}
MLAPILinearCombination_LDADD        = ${TEST_LIBS}
