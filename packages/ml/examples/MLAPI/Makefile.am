#@HEADER
#@HEADER

# include definitions of ML_EXPORT_INCLUDES and ML_EXPORT_LIBS
include $(top_builddir)/Makefile.export

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = PowerMethod          \
                  CreateTwoLevelDD     \
                  Preconditioner       \
		  Compare

PowerMethod_SOURCES = $(srcdir)/PowerMethod.cpp
PowerMethod_CPPFLAGS = ${ml_common_include_locations}
PowerMethod_CXXFLAGS = ${ml_common_include_locations}
PowerMethod_LDADD = ${ml_libs} 
PowerMethod_LINK = ${mllink}

CreateTwoLevelDD_SOURCES = $(srcdir)/CreateTwoLevelDD.cpp
CreateTwoLevelDD_CPPFLAGS = ${ml_common_include_locations}
CreateTwoLevelDD_CXXFLAGS = ${ml_common_include_locations}
CreateTwoLevelDD_LDADD = ${ml_libs} 
CreateTwoLevelDD_LINK = ${mllink}

Preconditioner_SOURCES = $(srcdir)/Preconditioner.cpp
Preconditioner_CPPFLAGS = ${ml_common_include_locations}
Preconditioner_CXXFLAGS = ${ml_common_include_locations}
Preconditioner_LDADD = ${ml_libs} 
Preconditioner_LINK = ${mllink}

Compare_SOURCES = $(srcdir)/Compare.cpp
Compare_CPPFLAGS = ${ml_common_include_locations}
Compare_CXXFLAGS = ${ml_common_include_locations}
Compare_LDADD = ${ml_libs} 
Compare_LINK = ${mllink}

# source file definitions, link lines
include ../Makefile-common.include
