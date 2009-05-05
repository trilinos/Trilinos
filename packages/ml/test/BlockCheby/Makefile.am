# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
include $(top_builddir)/Makefile.export.ml

TEST_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
TEST_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

EXTRA_DIST = samplemat.dat.gz blockids.dat.gz localids_in_blocks.dat.gz

noinst_PROGRAMS = BlockCheby GetNecessaryFiles

BlockCheby_SOURCES      = $(srcdir)/BlockCheby.cpp
BlockCheby_DEPENDENCIES = $(top_builddir)/src/libml.a
BlockCheby_CXXFLAGS     = ${TEST_INCLUDES}
BlockCheby_LDADD        = ${TEST_LIBS}


GetNecessaryFiles_SOURCES= $(top_srcdir)/test/BlockCheby

#copy these files over to the build area
GetNecessaryFiles$(EXEEXT): $(GetNecessaryFiles_SOURCES)
	@if test ! -e ./blockids.dat; then cp -f $(srcdir)/blockids.dat.gz .; gunzip -f ./blockids.dat.gz ; fi
	@if test ! -e ./localids_in_blocks.dat; then cp -f $(srcdir)/localids_in_blocks.dat.gz .; gunzip -f ./localids_in_blocks.dat.gz ; fi
	@if test ! -e ./samplemat.dat; then cp -f $(srcdir)/samplemat.dat.gz .; gunzip -f samplemat.dat.gz ; fi
