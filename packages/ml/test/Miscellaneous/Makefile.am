# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
include $(top_builddir)/Makefile.export.ml

TEST_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
TEST_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

noinst_PROGRAMS = ml_Miscellaneous

ml_Miscellaneous_SOURCES      = $(srcdir)/cxx_main.cpp
ml_Miscellaneous_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_Miscellaneous_CXXFLAGS     = ${TEST_INCLUDES}
ml_Miscellaneous_LDADD        = ${TEST_LIBS}
