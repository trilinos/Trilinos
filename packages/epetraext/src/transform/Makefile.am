# @HEADER
# @HEADER


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside EpetraExt. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

TRANSFORM_H = \
	$(srcdir)/Epetra_Transform.h \
	$(srcdir)/EDT_CrsGraph_View.h \
	$(srcdir)/EDT_CrsMatrix_View.h \
	$(srcdir)/EDT_MultiVector_View.h \
	$(srcdir)/EDT_CrsGraph_Overlap.h \
	$(srcdir)/EDT_CrsGraph_Transpose.h \
        $(srcdir)/EDT_CrsGraph_MapColoring.h \
        $(srcdir)/EDT_CrsGraph_MapColoringIndex.h \
	$(srcdir)/EDT_LinearProblem_StaticCondensation.h

TRANSFORM = \
	$(srcdir)/EDT_CrsGraph_View.cpp \
	$(srcdir)/EDT_CrsMatrix_View.cpp \
	$(srcdir)/EDT_MultiVector_View.cpp \
	$(srcdir)/EDT_CrsGraph_Overlap.cpp \
	$(srcdir)/EDT_CrsGraph_Transpose.cpp \
        $(srcdir)/EDT_CrsGraph_MapColoring.cpp \
        $(srcdir)/EDT_CrsGraph_MapColoringIndex.cpp \
	$(srcdir)/EDT_LinearProblem_StaticCondensation.cpp

BTF_H = \
        $(srddir)/ET_CrsGraph_BTF.h

BTF = \
        $(srddir)/ET_CrsGraph_BTF.cpp \
        $(srddir)/btf/concmp.f \
        $(srddir)/btf/genbtf.f \
        $(srddir)/btf/corsum.f \
        $(srddir)/btf/finsum.f \
        $(srddir)/btf/fnrsum.f \
        $(srddir)/btf/izero.f \
        $(srddir)/btf/maxmatch.f \
        $(srddir)/btf/rectblk.f \
        $(srddir)/btf/mmc13e.f \
        $(srddir)/btf/prtivs.f \
        $(srddir)/btf/strchk.f \
        $(srddir)/btf/mattrans.f

EXTRA_MPI_H = \
	$(srcdir)/EDT_CrsGraph_Zoltan.h \
	$(srcdir)/EDT_CrsGraph_ZoltanOrder.h \
	$(srcdir)/Epetra_ZoltanQuery.h

EXTRA_MPI = \
	$(srcdir)/EDT_CrsGraph_Zoltan.cpp \
	$(srcdir)/EDT_CrsGraph_ZoltanOrder.cpp \
	$(srcdir)/Epetra_ZoltanQuery.cpp

if HAVE_MPI

MPI_H = \
        $(EXTRA_MPI_H)
MPI = \
        $(EXTRA_MPI)

else

MPI_H =
MPI =

endif

# ------------------------------------------------------------------------
# EpetraExt library specifications
# ------------------------------------------------------------------------

EPETRAEXT_LIB = libepetraext.a

EPETRAEXT_H = \
	$(TRANSFORM_H) \
	$(MPI_H)

libepetraext_a_SOURCES = \
	$(TRANSFORM) \
	$(MPI)

#We need the following to support separate build directories
AM_CPPFLAGS=-I$(top_builddir)/../epetra/src -I$(top_srcdir)/../epetra/src \
            -I$(top_srcdir)/../../external/Zoltan_C++
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libepetraext_a_AR = $(ALTERNATE_AR)
else

libepetraext_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# On alpha-dec machines, the extra .o's are needed for Epetra only.
# ------------------------------------------------------------------------

libepetraext_a_LIBADD = $(EPETRAEXT_TEMPLATE_OBJS) $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(EPETRAEXT_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(pkgincludedir)
# ------------------------------------------------------------------------

pkginclude_HEADERS = $(EPETRAEXT_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(pkgexecincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: Epetra_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

pkgexecincludedir = $(exec_prefix)/include/epetraext
nodist_pkgexecinclude_HEADERS =../EpetraExt_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 



