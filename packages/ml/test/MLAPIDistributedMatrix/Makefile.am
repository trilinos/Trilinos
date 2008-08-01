# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
include $(top_builddir)/Makefile.export.ml

TEST_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
TEST_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

noinst_PROGRAMS = MLAPIDistributedMatrix

MLAPIDistributedMatrix_SOURCES      = $(srcdir)/cxx_main.cpp
MLAPIDistributedMatrix_DEPENDENCIES = $(top_builddir)/src/libml.a
MLAPIDistributedMatrix_CXXFLAGS     = ${TEST_INCLUDES}
MLAPIDistributedMatrix_LDADD        = ${TEST_LIBS}
