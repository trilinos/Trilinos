# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
include $(top_builddir)/Makefile.export.ml

TEST_LIBS     = $(ML_LIBS)
TEST_INCLUDES = $(ML_INCLUDES)

EXEEXT = .exe

noinst_PROGRAMS = AztecSimple

CCLD = $(CXX)

AztecSimple_SOURCES      = $(srcdir)/AztecSimple.c
AztecSimple_DEPENDENCIES = $(top_builddir)/src/libml.a
AztecSimple_CFLAGS       = ${TEST_INCLUDES}
AztecSimple_LDADD        = ${TEST_LIBS}
