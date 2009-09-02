# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
include $(top_builddir)/Makefile.export.ml

TEST_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
TEST_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES) $(ISORROPIA_INCLUDES)

EXEEXT = .exe

noinst_PROGRAMS = MultiLevelPreconditioner_Reduce

MultiLevelPreconditioner_Reduce_SOURCES      = $(srcdir)/MultiLevelPreconditioner_Reduce.cpp
MultiLevelPreconditioner_Reduce_DEPENDENCIES = $(top_builddir)/src/libml.a
MultiLevelPreconditioner_Reduce_CXXFLAGS     = ${TEST_INCLUDES}
MultiLevelPreconditioner_Reduce_LDADD        = ${TEST_LIBS}
