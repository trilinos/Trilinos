# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */

include $(top_builddir)/Makefile.export.ml

TEST_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
TEST_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

noinst_PROGRAMS = MLAPICompare

MLAPICompare_SOURCES      = $(srcdir)/cxx_main.cpp
MLAPICompare_DEPENDENCIES = $(top_builddir)/src/libml.a
MLAPICompare_CXXFLAGS     = ${TEST_INCLUDES}
MLAPICompare_LDADD        = ${TEST_LIBS}
