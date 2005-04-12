#@HEADER
#@HEADER

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export
include ../Makefile-common.include

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = Laplacian3D AdvDiff2D MLAPI TestAll

Laplacian3D_SOURCES = $(srcdir)/Laplacian3D.cpp
Laplacian3D_CPPFLAGS = ${ml_common_include_locations}
Laplacian3D_CXXFLAGS = ${ml_common_include_locations}
Laplacian3D_LDADD = ${ml_libs} 
Laplacian3D_LINK = ${mllink}
Laplacian3D_DEPENDENCIES = $(top_builddir)/src/libml.a

AdvDiff2D_SOURCES = $(srcdir)/AdvDiff2D.cpp
AdvDiff2D_CPPFLAGS = ${ml_common_include_locations}
AdvDiff2D_CXXFLAGS = ${ml_common_include_locations}
AdvDiff2D_LDADD = ${ml_libs} 
AdvDiff2D_LINK = ${mllink}
AdvDiff2D_DEPENDENCIES = $(top_builddir)/src/libml.a

MLAPI_SOURCES = $(srcdir)/MLAPI.cpp
MLAPI_CPPFLAGS = ${ml_common_include_locations}
MLAPI_CXXFLAGS = ${ml_common_include_locations}
MLAPI_LDADD = ${ml_libs} 
MLAPI_LINK = ${mllink}
MLAPI_DEPENDENCIES = $(top_builddir)/src/libml.a

TestAll_SOURCES = $(srcdir)/TestAll.cpp
TestAll_CPPFLAGS = ${ml_common_include_locations}
TestAll_CXXFLAGS = ${ml_common_include_locations}
TestAll_LDADD = ${ml_libs} 
TestAll_LINK = ${mllink}
TestAll_DEPENDENCIES = $(top_builddir)/src/libml.a
