# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
include $(top_builddir)/Makefile.export.ml

TEST_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
TEST_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

noinst_PROGRAMS = MultiLevelPreconditioner_Restart

MultiLevelPreconditioner_Restart_SOURCES      = $(srcdir)/MultiLevelPreconditioner_Restart.cpp
MultiLevelPreconditioner_Restart_DEPENDENCIES = $(top_builddir)/src/libml.a
MultiLevelPreconditioner_Restart_CXXFLAGS     = ${TEST_INCLUDES}
MultiLevelPreconditioner_Restart_LDADD        = ${TEST_LIBS}
