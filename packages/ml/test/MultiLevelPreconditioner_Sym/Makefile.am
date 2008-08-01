# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
include $(top_builddir)/Makefile.export.ml

TEST_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
TEST_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

noinst_PROGRAMS = MultiLevelPreconditioner_Sym 

MultiLevelPreconditioner_Sym_SOURCES      = $(srcdir)/MultiLevelPreconditioner_Sym.cpp
MultiLevelPreconditioner_Sym_DEPENDENCIES = $(top_builddir)/src/libml.a
MultiLevelPreconditioner_Sym_CXXFLAGS     = ${TEST_INCLUDES}
MultiLevelPreconditioner_Sym_LDADD        = ${TEST_LIBS}
