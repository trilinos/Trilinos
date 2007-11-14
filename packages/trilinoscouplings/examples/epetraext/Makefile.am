#@HEADER
#@HEADER

# include definitions of TRILINOSCOUPLINGS_INCLUDES and TRILINOSCOUPLINGS_LIBS
include $(top_builddir)/Makefile.export.trilinoscouplings

EXEEXT = .exe

# List the names of all exe's in this directory here

noinst_PROGRAMS = epetraext_isorropia_lptrans_example

epetraext_isorropia_lptrans_example_SOURCES      = $(srcdir)/EpetraExt_Isorropia_LPTrans_Ex.cpp 
epetraext_isorropia_lptrans_example_CXXFLAGS     = $(TRILINOSCOUPLINGS_INCLUDES)
epetraext_isorropia_lptrans_example_LDADD        = $(TRILINOSCOUPLINGS_LIBS)
epetraext_isorropia_lptrans_example_DEPENDENCIES = $(TRILINOSCOUPLINGS_LIBS)
