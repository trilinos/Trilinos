
# include definitions of TRILINOSCOUPLINGS_INCLUDES and TRILINOSCOUPLINGS_LIBS
include $(top_builddir)/Makefile.export.trilinoscouplings

EXEEXT = .exe

# List the names of all exe's in this directory here

if BUILD_TPETRA_ZOLTAN
ZOLTAN_EX = tpetra_zoltan_lptrans_example
else
ZOLTAN_EX = 
endif

noinst_PROGRAMS = $(ZOLTAN_EX)

tpetra_zoltan_lptrans_example_SOURCES      = $(srcdir)/EEP_Tpetra_Zoltan_LPTrans_Ex.cpp 
tpetra_zoltan_lptrans_example_CXXFLAGS     = $(TRILINOSCOUPLINGS_INCLUDES)
tpetra_zoltan_lptrans_example_LDADD        = $(TRILINOSCOUPLINGS_LIBS)
tpetra_zoltan_lptrans_example_DEPENDENCIES = $(TRILINOSCOUPLINGS_LIBS)
