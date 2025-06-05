
# include definitions of TRILINOSCOUPLINGS_INCLUDES and TRILINOSCOUPLINGS_LIBS
include $(top_builddir)/Makefile.export.trilinoscouplings

EXEEXT = .exe

# List the names of all exe's in this directory here

if BUILD_EPETRAEXT_AMESOS
AMESOS_EX = epetraext_amesosbtf_graphtrans_example epetraext_amesosbtf_matrixtrans_example
else
AMESOS_EX =
endif

noinst_PROGRAMS = $(AMESOS_EX)

epetraext_amesosbtf_graphtrans_example_SOURCES      = $(srcdir)/EpetraExt_AmesosBTF_GraphTrans_Ex.cpp 
epetraext_amesosbtf_graphtrans_example_CXXFLAGS     = $(TRILINOSCOUPLINGS_INCLUDES)
epetraext_amesosbtf_graphtrans_example_LDADD        = $(TRILINOSCOUPLINGS_LIBS)
epetraext_amesosbtf_graphtrans_example_DEPENDENCIES = $(TRILINOSCOUPLINGS_LIBS)

epetraext_amesosbtf_matrixtrans_example_SOURCES      = $(srcdir)/EpetraExt_AmesosBTF_MatrixTrans_Ex.cpp 
epetraext_amesosbtf_matrixtrans_example_CXXFLAGS     = $(TRILINOSCOUPLINGS_INCLUDES)
epetraext_amesosbtf_matrixtrans_example_LDADD        = $(TRILINOSCOUPLINGS_LIBS)
epetraext_amesosbtf_matrixtrans_example_DEPENDENCIES = $(TRILINOSCOUPLINGS_LIBS)
