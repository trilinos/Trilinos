# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
include $(top_builddir)/Makefile.export.ml

TEST_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
TEST_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXTRA_DIST= EvaluateDefinitions baseline1

GetNecessaryFiles_SOURCES= $(top_srcdir)/test/METISCoarsening
#copy these files over to the build area
GetNecessaryFiles$(EXEEXT): $(GetNecessaryFiles_SOURCES)
	@if test ! -e ./EvaluateDefinitions; then cp -f $(srcdir)/EvaluateDefinitions .; fi
	@if test ! -e ./baseline1 ; then cp -f $(srcdir)/baseline1 .; fi

EXEEXT = .exe

noinst_PROGRAMS = METISCoarsening GetNecessaryFiles

METISCoarsening_SOURCES      = $(srcdir)/cxx_main.cpp
METISCoarsening_DEPENDENCIES = $(top_builddir)/src/libml.a
METISCoarsening_CXXFLAGS     = $(TEST_INCLUDES)
METISCoarsening_LDADD        = ${TEST_LIBS}
