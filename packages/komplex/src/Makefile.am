# @HEADER
# @HEADER


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside Ifpack. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------
# For now, these are just guesses as to where the files belong.

CORE_H = \
	$(srcdir)/azk_komplex.h

CORE = \
	$(srcdir)/azk_create_linsys.c \
	$(srcdir)/azk_create_precon.c \
	$(srcdir)/azk_destroy_linsys.c \
	$(srcdir)/azk_destroy_precon.c \
	$(srcdir)/azk_create_matrix.c \
	$(srcdir)/azk_create_vector.c \
	$(srcdir)/azk_destroy_matrix.c \
	$(srcdir)/azk_destroy_vector.c 

### add in paths as needed - trying to prevent 
### unnecessary paths from being listed - see old.Makefile for ideas

AM_CPPFLAGS =-I$(top_srcdir)/../aztecoo/src -I$(top_builddir)/../aztecoo/src 
#	-I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src 

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

DIST_H = 

DIST = 

COMP_H = 

COMP = \
	$(srcdir)/azk_extract_solution.c \
	$(srcdir)/azk_matvec.c \
	$(srcdir)/azk_precon.c \
	$(srcdir)/azk_permute_ri.c \
	$(srcdir)/azk_residual_norm.c 



# ------------------------------------------------------------------------
# Komplex library specifications
# ------------------------------------------------------------------------

KOMPLEX_LIB = libkomplex.a


KOMPLEX_H = \
	$(CORE_H) \
	$(SERIAL_H) \
	$(USERBASE_H) \
	$(LANG_WRAPS_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(MPI_H) \
	$(MPISMP_H)

libkomplex_a_SOURCES = \
	$(CORE) \
	$(SERIAL) \
	$(LANG_WRAPS) \
	$(DIST) \
	$(COMP) \
	$(MPI) \
	$(MPISMP)

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libkomplex_a_AR = $(ALTERNATE_AR)
else

libkomplex_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libkomplex_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(KOMPLEX_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(pkgincludedir)
# ------------------------------------------------------------------------

pkginclude_HEADERS = $(KOMPLEX_H)

# ------------------------------------------------------------------------
# Special stuff to install in our special $(pkgexecincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: Ifpack_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

pkgexecincludedir = $(exec_prefix)/include
nodist_pkgexecinclude_HEADERS = Komplex_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






