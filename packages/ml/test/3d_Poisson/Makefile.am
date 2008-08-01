# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
EXEEXT = .exe

# include definitions of ML_INCLUDES and ML_LIBS
include $(top_builddir)/Makefile.export.ml

EXTRA_DIST= 3dPoissonEvaluateDefinitions baseline1 baseline2 baseline3 baseline5 baseline6 baseline7 ml_inputfile1 ml_inputfile2 ml_inputfile3 ml_inputfile5 ml_inputfile6 ml_inputfile7

noinst_PROGRAMS = 3d_Poisson GetNecessaryFiles

3d_Poisson_SOURCES = $(top_srcdir)/examples/Other/ml_example3d.c

3d_Poisson$(EXT): $(3d_Poisson_SOURCES)

GetNecessaryFiles_SOURCES= $(top_srcdir)/test/3d_Poisson
#copy these files over to the build area
GetNecessaryFiles$(EXEEXT): $(GetNecessaryFiles_SOURCES)
	@if test ! -e ./3dPoissonEvaluateDefinitions; then cp -f $(srcdir)/3dPoissonEvaluateDefinitions .; fi
	@if test ! -e ./baseline1 ; then cp -f $(srcdir)/baseline1 .; fi
	@if test ! -e ./baseline2 ; then cp -f $(srcdir)/baseline2 .; fi
	@if test ! -e ./baseline3 ; then cp -f $(srcdir)/baseline3 .; fi
	@if test ! -e ./baseline5 ; then cp -f $(srcdir)/baseline5 .; fi
	@if test ! -e ./baseline6 ; then cp -f $(srcdir)/baseline6 .; fi
	@if test ! -e ./baseline7 ; then cp -f $(srcdir)/baseline7 .; fi
	@if test ! -e ./ml_inputfile1 ; then cp -f $(srcdir)/ml_inputfile1 .; fi
	@if test ! -e ./ml_inputfile2 ; then cp -f $(srcdir)/ml_inputfile2 .; fi
	@if test ! -e ./ml_inputfile3 ; then cp -f $(srcdir)/ml_inputfile3 .; fi
	@if test ! -e ./ml_inputfile5 ; then cp -f $(srcdir)/ml_inputfile5 .; fi
	@if test ! -e ./ml_inputfile6 ; then cp -f $(srcdir)/ml_inputfile6 .; fi
	@if test ! -e ./ml_inputfile7 ; then cp -f $(srcdir)/ml_inputfile7 .; fi

3d_Poisson_DEPENDENCIES = \
    $(top_builddir)/src/libml.a

ml_common_include_locations = $(ML_INCLUDES)

ml_libs = $(ML_LIBS) \
          $(LIBS)           \
          $(LAPACK_LIBS)    \
          $(BLAS_LIBS)      \
          $(FLIBS)

mllink = $(CXX) $(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

3d_Poisson_CPPFLAGS = ${ml_common_include_locations} -DML_BENCHMARK
3d_Poisson_CXXFLAGS = ${ml_common_include_locations} -DML_BENCHMARK
3d_Poisson_LDADD = ${ml_libs}
3d_Poisson_LINK = ${mllink}
