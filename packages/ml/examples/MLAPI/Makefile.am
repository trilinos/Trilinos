#@HEADER
#@HEADER

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export
ml_common_include_locations = $(ML_EXPORT_INCLUDES)
ml_libs = $(ML_EXPORT_LIBS) $(LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
mllink = $(CXX) $(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = Blackboard           \
                  PowerMethod          \
                  TwoLevelDDAdditive   \
                  MultiLevel         

Blackboard_SOURCES  = $(srcdir)/Blackboard.cpp
Blackboard_CPPFLAGS = ${ml_common_include_locations}
Blackboard_CXXFLAGS = ${ml_common_include_locations}
Blackboard_LDADD    = ${ml_libs} 
Blackboard_LINK     = ${mllink}

PowerMethod_SOURCES  = $(srcdir)/PowerMethod.cpp
PowerMethod_CPPFLAGS = ${ml_common_include_locations}
PowerMethod_CXXFLAGS = ${ml_common_include_locations}
PowerMethod_LDADD    = ${ml_libs} 
PowerMethod_LINK     = ${mllink}

TwoLevelDDAdditive_SOURCES  = $(srcdir)/TwoLevelDDAdditive.cpp
TwoLevelDDAdditive_CPPFLAGS = ${ml_common_include_locations}
TwoLevelDDAdditive_CXXFLAGS = ${ml_common_include_locations}
TwoLevelDDAdditive_LDADD    = ${ml_libs} 
TwoLevelDDAdditive_LINK     = ${mllink}

MultiLevel_SOURCES  = $(srcdir)/MultiLevel.cpp
MultiLevel_CPPFLAGS = ${ml_common_include_locations}
MultiLevel_CXXFLAGS = ${ml_common_include_locations}
MultiLevel_LDADD    = ${ml_libs} 
MultiLevel_LINK     = ${mllink}
