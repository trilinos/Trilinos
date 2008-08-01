# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
include $(top_builddir)/Makefile.export.ml

TEST_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
TEST_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

if ML_USING_AZTECOO
if ML_USING_TEUCHOS
if ML_USING_GALERI
if ML_USING_EPETRAEXT
if ML_USING_EPETRA
noinst_PROGRAMS = MatrixFree 
else
noinst_PROGRAMS = 
endif
else
noinst_PROGRAMS = 
endif
else
noinst_PROGRAMS = 
endif
else
noinst_PROGRAMS = 
endif
else
noinst_PROGRAMS = 
endif

MatrixFree_SOURCES      = $(srcdir)/MatrixFree.cpp
MatrixFree_DEPENDENCIES = $(top_builddir)/src/libml.a
MatrixFree_CXXFLAGS     = ${TEST_INCLUDES}
MatrixFree_LDADD        = ${TEST_LIBS}
