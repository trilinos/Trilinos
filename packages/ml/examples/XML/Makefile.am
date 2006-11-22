#@HEADER

include $(top_builddir)/Makefile.export.ml

EXAMPLE_LIBS     = $(ML_EXTRA_LIBS)     $(ML_LIBS)
EXAMPLE_INCLUDES = $(ML_EXTRA_INCLUDES) $(ML_INCLUDES)

EXEEXT = .exe

EXTRA_DIST = ml_ParameterList.xml

# List the names of all exe's in this directory here

noinst_PROGRAMS = ml_XML 

CCLD = $(CXX)

ml_XML_SOURCES      = $(srcdir)/ml_XML.cpp
ml_XML_DEPENDENCIES = $(top_builddir)/src/libml.a
ml_XML_CXXFLAGS     = $(EXAMPLE_INCLUDES)
ml_XML_LDADD        = $(EXAMPLE_LIBS)
