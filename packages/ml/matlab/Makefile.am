# ************************************************************************* */
# See the file COPYRIGHT for a complete copyright notice, contact person,   */
# and disclaimer.                                                           */
# ************************************************************************* */
# Include export Makefiles
include  $(top_builddir)/Makefile.export.ml

# I need extra_libs because of MLAPI::Gallery, which calls Galeri.
ML_MLMEX_INCLUDES = $(ML_INCLUDES) $(ML_EXTRA_INCLUDES) $(ML_MATLAB_INCLUDES)
ML_MLMEX_LIBS     = $(ML_LIBS) $(ML_EXTRA_LIBS)

# The MEXEXT_VALUE is set in AZ_MATLAB_ROOT from the mexext script, the only portable
# way of figuring out what my mex system extension is.  Of course, this means cross-compiling
# probably won't work.  C'est la vie.
EXEEXT = .$(MEXEXT_VALUE)
EXTRA_DIST = ml.m

if BUILD_ML_MATLAB
bin_PROGRAMS = mlmex 
bin_SCRIPTS = ml.m

mlmex_SOURCES      = mlmex.cpp
mlmex_DEPENDENCIES = $(top_builddir)/src/libml.a
mlmex_CXXFLAGS     = ${ML_MLMEX_INCLUDES} -fPIC
mlmex_LDADD        = ${ML_MLMEX_LIBS} 
mlmex_LINK         = $(MEX) -cxx -o $@
endif



