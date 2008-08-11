# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
include $(top_builddir)/Makefile.export.ml

TEST_LIBS = $(ML_LIBS)  $(ML_EXTRA_LIBS) 
TEST_INCLUDES = $(ML_INCLUDES) $(ML_EXTRA_INCLUDES)

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_periodic

ml_periodic_SOURCES      = $(srcdir)/ml_periodic.cpp
ml_periodic_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_periodic_CXXFLAGS     = $(TEST_INCLUDES)
ml_periodic_LDADD        = $(TEST_LIBS)
