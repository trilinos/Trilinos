# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
include $(top_builddir)/Makefile.export.ml

TEST_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
TEST_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

noinst_PROGRAMS = MLAPISerialMatrix

MLAPISerialMatrix_SOURCES      = $(srcdir)/cxx_main.cpp
MLAPISerialMatrix_DEPENDENCIES = $(top_builddir)/src/libml.a
MLAPISerialMatrix_CXXFLAGS     = ${TEST_INCLUDES}
MLAPISerialMatrix_LDADD        = ${TEST_LIBS}
