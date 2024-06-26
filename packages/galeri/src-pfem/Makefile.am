
# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package galeri subdirectory src-pfem built successfully."
	@echo ""

CORE_H = $(srcdir)/Galeri_core_Object.h \
	 $(srcdir)/Galeri_grid_Triangle.h \
	 $(srcdir)/Galeri_core_Workspace.h \
	 $(srcdir)/Galeri_problem_Base.h \
	 $(srcdir)/Galeri_grid_Element.h \
	 $(srcdir)/Galeri_problem_ScalarLaplacian.h \
	 $(srcdir)/Galeri_grid_Generator.h \
	 $(srcdir)/Galeri_problem_VectorLaplacian.h \
	 $(srcdir)/Galeri_grid_Hex.h \
	 $(srcdir)/Galeri_quadrature_Element.h \
	 $(srcdir)/Galeri_grid_Loadable.h \
	 $(srcdir)/Galeri_quadrature_Hex.h \
	 $(srcdir)/Galeri_grid_Point.h \
	 $(srcdir)/Galeri_quadrature_Quad.h \
	 $(srcdir)/Galeri_grid_Quad.h \
	 $(srcdir)/Galeri_quadrature_Segment.h \
	 $(srcdir)/Galeri_grid_Rebalance.h \
	 $(srcdir)/Galeri_quadrature_Tet.h \
	 $(srcdir)/Galeri_grid_Segment.h \
	 $(srcdir)/Galeri_quadrature_Triangle.h \
	 $(srcdir)/Galeri_grid_SerialXML.h \
	 $(srcdir)/Galeri_viz_MEDIT.h \
	 $(srcdir)/Galeri_grid_Tet.h \
	 $(srcdir)/Galeri_viz_VTK.h

CORE =  $(srcdir)/Galeri_core_Workspace.cpp \
	$(srcdir)/Galeri_grid_Generator.cpp \
	$(srcdir)/Galeri_grid_Loadable.cpp

include $(top_builddir)/Makefile.export.galeri

AM_CPPFLAGS = $(GALERI_INCLUDES)

# ------------------------------------------------------------------------
# pfem library specifications
# ------------------------------------------------------------------------
PFEM_LIB = libpfem.a

PFEM_H = $(CORE_H) 

libpfem_a_SOURCES = $(CORE) 

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR
libpfem_a_AR = $(ALTERNATE_AR)
else
libpfem_a_AR = $(AR) cru
endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libpfem_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(PFEM_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(PFEM_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# All Trilinos headers are now installed in the same directory
execincludedir = $(exec_prefix)/include

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 
