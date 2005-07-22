#@HEADER

# include definitions of ML_INCLUDES and ML_LIBS
include $(top_builddir)/Makefile.export.ml

if USING_GNUMAKE
EXAMPLE_LIBS     = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(ML_LIBS) $(ML_EXTRA_LIBS))
EXAMPLE_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(ML_INCLUDES) $(ML_EXTRA_INCLUDES))
else
EXAMPLE_LIBS     = $(ML_LIBS)     $(ML_EXTRA_LIBS)
EXAMPLE_INCLUDES = $(ML_INCLUDES) $(ML_EXTRA_INCLUDES)
endif

AM_CPPFLAGS = $(EXAMPLE_INCLUDES)

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = Blackboard            \
                  PowerMethod           \
                  TwoLevelDDAdditive    \
                  MultiLevel            \
                  Richardson            \
                  AdaptiveSA            \
                  EpetraInterface

EXTRA_DIST = README

Blackboard_SOURCES      = $(srcdir)/Blackboard.cpp
Blackboard_DEPENDENCIES = $(top_builddir)/src/libml.a
Blackboard_LDADD        = ${EXAMPLE_LIBS} 

PowerMethod_SOURCES      = $(srcdir)/PowerMethod.cpp
PowerMethod_DEPENDENCIES = $(top_builddir)/src/libml.a
PowerMethod_LDADD        = ${EXAMPLE_LIBS} 

TwoLevelDDAdditive_SOURCES      = $(srcdir)/TwoLevelDDAdditive.cpp
TwoLevelDDAdditive_DEPENDENCIES = $(top_builddir)/src/libml.a
TwoLevelDDAdditive_LDADD        = ${EXAMPLE_LIBS} 

MultiLevel_SOURCES      = $(srcdir)/MultiLevel.cpp
MultiLevel_DEPENDENCIES = $(top_builddir)/src/libml.a
MultiLevel_LDADD        = ${EXAMPLE_LIBS} 

AdaptiveSA_SOURCES      = $(srcdir)/AdaptiveSA.cpp
AdaptiveSA_DEPENDENCIES = $(top_builddir)/src/libml.a
AdaptiveSA_LDADD        = ${EXAMPLE_LIBS} 

EpetraInterface_SOURCES      = $(srcdir)/EpetraInterface.cpp
EpetraInterface_DEPENDENCIES = $(top_builddir)/src/libml.a
EpetraInterface_LDADD        = ${EXAMPLE_LIBS} 

Richardson_SOURCES      = $(srcdir)/Richardson.cpp
Richardson_DEPENDENCIES = $(top_builddir)/src/libml.a
Richardson_LDADD        = ${EXAMPLE_LIBS} 
