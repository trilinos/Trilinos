# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
include $(top_builddir)/Makefile.export.ml

TEST_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
TEST_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

noinst_PROGRAMS = ML_Operator2Epetra_RowMatrix 

ML_Operator2Epetra_RowMatrix_SOURCES      = $(srcdir)/ML_Operator2Epetra_RowMatrix.cpp

ML_Operator2Epetra_RowMatrix_DEPENDENCIES = $(top_builddir)/src/libml.a
ML_Operator2Epetra_RowMatrix_CXXFLAGS     = ${TEST_INCLUDES}
ML_Operator2Epetra_RowMatrix_LDADD        = ${TEST_LIBS}
