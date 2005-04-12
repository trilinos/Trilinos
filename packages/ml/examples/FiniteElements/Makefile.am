#@HEADER
#@HEADER

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export
include ../Makefile-common.include

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

Laplacian3D_SOURCES = $(srcdir)/Laplacian3D.cpp $(FEHEADERS)
Laplacian3D_CPPFLAGS = ${ml_common_include_locations}
Laplacian3D_CXXFLAGS = ${ml_common_include_locations}
Laplacian3D_LDADD = ${ml_libs} 
Laplacian3D_LINK = ${mllink}
Laplacian3D_DEPENDENCIES = $(top_builddir)/src/libml.a

AdvDiff2D_SOURCES = $(srcdir)/AdvDiff2D.cpp $(FEHEADERS)
AdvDiff2D_CPPFLAGS = ${ml_common_include_locations}
AdvDiff2D_CXXFLAGS = ${ml_common_include_locations}
AdvDiff2D_LDADD = ${ml_libs} 
AdvDiff2D_LINK = ${mllink}
AdvDiff2D_DEPENDENCIES = $(top_builddir)/src/libml.a

MLAPI_SOURCES = $(srcdir)/MLAPI.cpp $(FEHEADERS)
MLAPI_CPPFLAGS = ${ml_common_include_locations}
MLAPI_CXXFLAGS = ${ml_common_include_locations}
MLAPI_LDADD = ${ml_libs} 
MLAPI_LINK = ${mllink}
MLAPI_DEPENDENCIES = $(top_builddir)/src/libml.a

TestAll_SOURCES = $(srcdir)/TestAll.cpp $(FEHEADERS)
TestAll_CPPFLAGS = ${ml_common_include_locations}
TestAll_CXXFLAGS = ${ml_common_include_locations}
TestAll_LDADD = ${ml_libs} 
TestAll_LINK = ${mllink}
TestAll_DEPENDENCIES = $(top_builddir)/src/libml.a
