#@HEADER
#@HEADER

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export
include ../Makefile-common.include

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = Blackboard            \
                  PowerMethod           \
                  TwoLevelDDAdditive    \
                  MultiLevel            \
                  Richardson            \
                  AdaptiveSA            \
                  AdaptiveSA_Elasticity \
                  AdaptiveSA_Tire       \
                  AdaptiveSA_SAMIS      \
                  NonSymmetricSA

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

AdaptiveSA_SOURCES  = $(srcdir)/AdaptiveSA.cpp
AdaptiveSA_CPPFLAGS = ${ml_common_include_locations}
AdaptiveSA_CXXFLAGS = ${ml_common_include_locations}
AdaptiveSA_LDADD    = ${ml_libs} 
AdaptiveSA_LINK     = ${mllink}

AdaptiveSA_Elasticity_SOURCES  = $(srcdir)/AdaptiveSA_Elasticity.cpp
AdaptiveSA_Elasticity_CPPFLAGS = ${ml_common_include_locations}
AdaptiveSA_Elasticity_CXXFLAGS = ${ml_common_include_locations}
AdaptiveSA_Elasticity_LDADD    = ${ml_libs} 
AdaptiveSA_Elasticity_LINK     = ${mllink}

AdaptiveSA_Tire_SOURCES  = $(srcdir)/AdaptiveSA_Tire.cpp
AdaptiveSA_Tire_CPPFLAGS = ${ml_common_include_locations}
AdaptiveSA_Tire_CXXFLAGS = ${ml_common_include_locations}
AdaptiveSA_Tire_LDADD    = ${ml_libs} 
AdaptiveSA_Tire_LINK     = ${mllink}

Richardson_SOURCES  = $(srcdir)/Richardson.cpp
Richardson_CPPFLAGS = ${ml_common_include_locations}
Richardson_CXXFLAGS = ${ml_common_include_locations}
Richardson_LDADD    = ${ml_libs} 
Richardson_LINK     = ${mllink}

NonSymmetricSA_SOURCES  = $(srcdir)/NonSymmetricSA.cpp
NonSymmetricSA_CPPFLAGS = ${ml_common_include_locations}
NonSymmetricSA_CXXFLAGS = ${ml_common_include_locations}
NonSymmetricSA_LDADD    = ${ml_libs} 
NonSymmetricSA_LINK     = ${mllink}

AdaptiveSA_SAMIS_SOURCES  = $(srcdir)/AdaptiveSA_SAMIS.cpp
AdaptiveSA_SAMIS_CPPFLAGS = ${ml_common_include_locations}
AdaptiveSA_SAMIS_CXXFLAGS = ${ml_common_include_locations}
AdaptiveSA_SAMIS_LDADD    = ${ml_libs} 
AdaptiveSA_SAMIS_LINK     = ${mllink}
