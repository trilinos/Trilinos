#@HEADER
#@HEADER

# include definitions of ML_INCLUDES and ML_LIBS
include $(top_builddir)/Makefile.export.ml

if USING_GNUMAKE
EXAMPLE_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(ML_LIBS))
EXAMPLE_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(ML_INCLUDES))
else
EXAMPLE_LIBS     = $(ML_LIBS)
EXAMPLE_INCLUDES = $(ML_INCLUDES)
endif

AM_CPPFLAGS = $(EXAMPLE_INCLUDES)

EXEEXT = .exe

# List the names of all exe's in this directory here

FEHEADERS = \
ml_AbstractGrid.h               ml_QuadQuadrature.h        \
ml_AbstractProblem.h            ml_QuadRectangleGrid.h     \
ml_AbstractQuadrature.h         ml_SUPGVariational.h       \
ml_AbstractVariational.h        ml_TRIANGLEGrid.h          \
ml_BasicGrid.h                  ml_TetCubeGrid.h           \
ml_GalerkinVariational.h        ml_TetQuadrature.h         \
ml_HexCubeGrid.h                ml_TriangleQuadrature.h    \
ml_HexQuadrature.h              ml_TriangleRectangleGrid.h \
ml_LinearProblem.h              ml_Utils.h                 \
ml_MEDITInterface.h             ml_Workspace.h

noinst_PROGRAMS = Laplacian3D AdvDiff2D MLAPI TestAll

Laplacian3D_SOURCES      = $(srcdir)/Laplacian3D.cpp $(FEHEADERS)
Laplacian3D_DEPENDENCIES = $(top_builddir)/src/libml.a
Laplacian3D_LDADD        = ${EXAMPLE_LIBS} 

AdvDiff2D_SOURCES      = $(srcdir)/AdvDiff2D.cpp $(FEHEADERS)
AdvDiff2D_DEPENDENCIES = $(top_builddir)/src/libml.a
AdvDiff2D_LDADD        = ${EXAMPLE_LIBS} 

MLAPI_SOURCES      = $(srcdir)/MLAPI.cpp $(FEHEADERS)
MLAPI_DEPENDENCIES = $(top_builddir)/src/libml.a
MLAPI_LDADD        = ${EXAMPLE_LIBS} 

TestAll_SOURCES      = $(srcdir)/TestAll.cpp $(FEHEADERS)
TestAll_DEPENDENCIES = $(top_builddir)/src/libml.a
TestAll_LDADD        = ${EXAMPLE_LIBS} 
