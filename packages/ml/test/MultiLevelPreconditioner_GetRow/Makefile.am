# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
include $(top_builddir)/Makefile.export.ml

TEST_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
TEST_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

noinst_PROGRAMS = MultiLevelPreconditioner_GetRow

MultiLevelPreconditioner_GetRow_SOURCES      = $(srcdir)/MultiLevelPreconditioner_GetRow.cpp
MultiLevelPreconditioner_GetRow_DEPENDENCIES = $(top_builddir)/src/libml.a
MultiLevelPreconditioner_GetRow_CXXFLAGS     = ${TEST_INCLUDES}
MultiLevelPreconditioner_GetRow_LDADD        = ${TEST_LIBS}
