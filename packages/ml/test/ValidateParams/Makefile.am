# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
include $(top_builddir)/Makefile.export.ml

TEST_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
TEST_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

if ML_USING_TEUCHOS
if ML_USING_AZTECOO
if ML_USING_EPETRA
if ML_USING_GALERI
if ML_USING_AMESOS

noinst_PROGRAMS = ml_Validate GetNecessaryFiles

GetNecessaryFiles_SOURCES= $(top_srcdir)/test/ValidateParams

#copy these files over to the build area
GetNecessaryFiles$(EXEEXT): $(GetNecessaryFiles_SOURCES)
	@if test ! -e ./input1.xml; then cp -f $(srcdir)/input1.xml .; fi
	@if test ! -e ./baseline1.txt; then cp -f $(srcdir)/baseline1.txt .; fi
	@if test ! -e ./evaluateCriteria; then cp -f $(srcdir)/evaluateCriteria .; fi

ml_Validate_SOURCES      = $(srcdir)/cxx_main.cpp
ml_Validate_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_Validate_CXXFLAGS     = ${TEST_INCLUDES}
ml_Validate_LDADD        = ${TEST_LIBS}

endif
endif
endif
endif
endif
