# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
include $(top_builddir)/Makefile.export.ml

TEST_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
TEST_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

noinst_PROGRAMS = SelfSmoother

SelfSmoother_SOURCES      = $(srcdir)/cxx_main.cpp
SelfSmoother_DEPENDENCIES = $(top_builddir)/src/libml.a
SelfSmoother_CXXFLAGS     = $(TEST_INCLUDES)
SelfSmoother_LDADD        = ${TEST_LIBS}
