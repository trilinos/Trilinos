# @HEADER
# @HEADER


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside Ifpack. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------
# For now, these are just guesses as to where the files belong.
# Files that prove to be troublesome have been commented out and 
# placed at the bottom of the appropriate section

#Ifpack_CrsGraph.h
#Ifpack_CrsRick.h
#ifp_blas1.h
#
#
#

CORE_H = \
	$(srcdir)/Ifpack_ConfigDefs.h

CORE = 

### add in paths as needed - trying to prevent 
### unnecessary paths from being listed - see old.Makefile for ideas
### Should I add a dependence on AztecOO??  Library doesn't necessarily have
### to be built - just need some headers???  Ask Mike about headers that 
### may be needed outside of the various packages, vs those that will not
### be needed and if I should split them up.

AM_CPPFLAGS =-I$(top_srcdir)/../aztecoo/src -I$(top_builddir)/../aztecoo/src \
	-I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src 

SERIAL_H = 

SERIAL = 

if HAVE_MPI

MPI_H = 

MPI = 

else

MPI_H =
MPI =

endif

#Ifpack does not yet support threads
#if HAVE_THREADS

MPISMP_H = 

MPISMP = 

#else

#MPISMP_H =
#MPISMP =

#endif

USERBASE_H = 

LANG_WRAPS_H = 

LANG_WRAPS = 

EXTRA_DIST_H = \
	$(srcdir)/az_ifpack.h \
	$(srcdir)/ifp_BlockMat.h \
	$(srcdir)/ifp_BlockVec.h \
	$(srcdir)/ifp_DenseMat.h \
	$(srcdir)/ifp_SparseUtil.h \
	$(srcdir)/ifp_biluk.h \
	$(srcdir)/ifp_brelax.h \
	$(srcdir)/ifp_c_wrappers.h \
	$(srcdir)/ifp_Matrix.h \
	$(srcdir)/ifp_LocalMat.h \
	$(srcdir)/ifp_LocalPrecon.h \
	$(srcdir)/ifp_ifpack.h \
	$(srcdir)/ifp_blas3.h \
	$(srcdir)/ifp_lapackd.h \
	$(srcdir)/ifp_arch.h \
	$(srcdir)/ifp_GlobalPrecon.h \
	$(srcdir)/ifp_Precon.h \
	$(srcdir)/ifp_spblas.h

EXTRA_DIST = \
	$(srcdir)/az_ifpack_prec_create.c \
	$(srcdir)/az_ifpack_precon.c  \
	$(srcdir)/az_ifpack_iterate.c \
	$(srcdir)/az_ifpack_prec_destroy.c \
	$(srcdir)/az_ifpack_solve.c \
	$(srcdir)/ifp_BlockMat.cpp \
	$(srcdir)/ifp_BlockVec.cpp \
	$(srcdir)/ifp_DenseMat.cpp \
	$(srcdir)/ifp_SparseUtil.cpp \
	$(srcdir)/ifp_biluk.cpp \
	$(srcdir)/ifp_brelax.cpp \
	$(srcdir)/ifp_c_wrappers.cpp \
	$(srcdir)/ifp_spmm.cpp \
	$(srcdir)/ifp_spsm.cpp 

if BUILD_OLDIFPACK

DIST_H = $(EXTRA_DIST_H) 

DIST = $(EXTRA_DIST)

else

DIST_H = 
DIST = 

endif

COMP_H = \
	$(srcdir)/Ifpack_IlukGraph.h \
	$(srcdir)/Ifpack_CrsRiluk.h \
	$(srcdir)/Ifpack_ScalingType.h

COMP = \
	$(srcdir)/Ifpack_IlukGraph.cpp \
	$(srcdir)/Ifpack_CrsRiluk.cpp 


# ------------------------------------------------------------------------
# Ifpack library specifications
# ------------------------------------------------------------------------

IFPACK_LIB = libifpack.a


IFPACK_H = \
	$(CORE_H) \
	$(SERIAL_H) \
	$(USERBASE_H) \
	$(LANG_WRAPS_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(MPI_H) \
	$(MPISMP_H)

libifpack_a_SOURCES = \
	$(CORE) \
	$(SERIAL) \
	$(LANG_WRAPS) \
	$(DIST) \
	$(COMP) \
	$(MPI) \
	$(MPISMP)

EXTRA_libifpack_a_SOURCES = \
	$(EXTRA_DIST) \
	$(EXTRA_DIST_H) 

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libifpack_a_AR = $(ALTERNATE_AR)
else

libifpack_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libifpack_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(IFPACK_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(pkgincludedir)
# ------------------------------------------------------------------------

pkginclude_HEADERS = $(IFPACK_H)

# ------------------------------------------------------------------------
# Special stuff to install in our special $(pkgexecincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: Ifpack_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

pkgexecincludedir = $(exec_prefix)/include/ifpack
nodist_pkgexecinclude_HEADERS = Ifpack_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






