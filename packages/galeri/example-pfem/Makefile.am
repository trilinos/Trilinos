
all-local:
	@echo ""
	@echo "Trilinos package galeri subdirectory example-pfem built successfully."
	@echo ""

include $(top_builddir)/Makefile.export.galeri

EXEEXT = .exe

noinst_PROGRAMS = ex1 ex2 ex3

ex1_SOURCES      = $(srcdir)/ex1.cpp
ex1_DEPENDENCIES = $(top_builddir)/src-pfem/libpfem.a
ex1_CXXFLAGS     = $(GALERI_INCLUDES)
ex1_LDADD        = $(GALERI_LIBS) 

ex2_SOURCES      = $(srcdir)/ex2.cpp
ex2_DEPENDENCIES = $(top_builddir)/src-pfem/libpfem.a
ex2_CXXFLAGS     = $(GALERI_INCLUDES)
ex2_LDADD        = $(GALERI_LIBS) 

ex3_SOURCES      = $(srcdir)/ex3.cpp
ex3_DEPENDENCIES = $(top_builddir)/src-pfem/libpfem.a
ex3_CXXFLAGS     = $(GALERI_INCLUDES)
ex3_LDADD        = $(GALERI_LIBS) 
