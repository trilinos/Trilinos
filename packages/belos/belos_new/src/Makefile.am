# @HEADER
# @HEADER


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside Ifpack. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

CORE_H = \
	$(srcdir)/BelosBlockCG.hpp \
	$(srcdir)/BelosBlockGmres.hpp \
	$(srcdir)/BelosConfigDefs.hpp \
	$(srcdir)/BelosEpetraOperator.hpp \
	$(srcdir)/BelosIterativeSolver.hpp \
	$(srcdir)/BelosLinearProblemManager.hpp \
	$(srcdir)/BelosMultiVec.hpp \
	$(srcdir)/BelosOperator.hpp \
	$(srcdir)/BelosPetraInterface.hpp \
	$(srcdir)/BelosStatusTest.hpp \
	$(srcdir)/BelosStatusTestMaxIters.hpp \
	$(srcdir)/BelosStatusTestMaxRestarts.hpp \
	$(srcdir)/BelosStatusTestResNorm.hpp \
	$(srcdir)/BelosStatusType.hpp \
	$(srcdir)/BelosReturnType.hpp

CORE = 

UTIL = \
	$(srcdir)/BelosVersion.cpp

### add in paths as needed

AM_CPPFLAGS =-I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src
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

#if HAVE_THREADS

MPISMP_H = 

MPISMP = 

#else

#MPISMP_H =
#MPISMP =

#endif

LANG_WRAPS_H = 

LANG_WRAPS = 

DIST_H = 

DIST = 

COMP_H = \
	$(srcdir)/BelosBlockCG.hpp \
	$(srcdir)/BelosBlockGmres.hpp

COMP = 



# ------------------------------------------------------------------------
# Belos library specifications
# ------------------------------------------------------------------------

BELOS_LIB = libbelos.a


BELOS_H = \
	$(CORE_H) \
	$(SERIAL_H) \
	$(LANG_WRAPS_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(MPI_H) \
	$(MPISMP_H)

libbelos_a_SOURCES = \
	$(CORE) \
	$(UTIL) \
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

libbelos_a_AR = $(ALTERNATE_AR)
else

libbelos_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libbelos_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(BELOS_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(BELOS_H)

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: Belos_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = Belos_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 



