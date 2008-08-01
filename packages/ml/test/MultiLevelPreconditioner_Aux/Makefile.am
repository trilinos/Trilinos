# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
include $(top_builddir)/Makefile.export.ml

TEST_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
TEST_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

noinst_PROGRAMS = MultiLevelPreconditioner_Aux

MultiLevelPreconditioner_Aux_SOURCES      = $(srcdir)/MultiLevelPreconditioner_Aux.cpp
MultiLevelPreconditioner_Aux_DEPENDENCIES = $(top_builddir)/src/libml.a
MultiLevelPreconditioner_Aux_CXXFLAGS     = ${TEST_INCLUDES}
MultiLevelPreconditioner_Aux_LDADD        = ${TEST_LIBS}
