# @HEADER
# @HEADER


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside the Belos tests. 
# The first is the list of all source and any other header files.
# ------------------------------------------------------------------------

CORE_H = 

CORE = 

UTIL_H = \
	$(srcdir)/Util.h 

UTIL = \
	$(srcdir)/read_hb_matrix.cpp \
	$(srcdir)/distrib_msr_matrix.cpp 

SERIAL_H = 

SERIAL = 

EXTRA_MPI_H = 

EXTRA_MPI = 

#if HAVE_MPI

#MPI_H = \
#	$(EXTRA_MPI_H)
#MPI = \
#	$(EXTRA_MPI)

#else

MPI_H =
MPI =

#endif

EXTRA_MPISMP_H = 

EXTRA_MPISMP = 

#if HAVE_THREADS

#MPISMP_H = \
#	$(EXTRA_MPISMP_H)
#MPISMP = \
#	$(EXTRA_MPISMP)

#else

MPISMP_H =
MPISMP =

#endif

DIST_H = 

DIST = 

COMP_H = 

COMP = 

# ------------------------------------------------------------------------
# Belos Test library specifications
# ------------------------------------------------------------------------

BELOSUTILS_LIB = libbelosutils.a


BELOSUTILS_H = \
	$(CORE_H) \
	$(UTIL_H) \
	$(SERIAL_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(MPI_H) \
	$(MPISMP_H)

libbelosutils_a_SOURCES = \
	$(CORE) \
	$(UTIL) \
	$(SERIAL) \
	$(DIST) \
	$(COMP) \
	$(MPI) \
	$(MPISMP)

#We need the following to support separate build directories
AM_CPPFLAGS=-I$(top_srcdir)/src -I$(top_builddir)/src \
		-I$(top_srcdir)/../epetra/src \
		-I$(top_builddir)/../epetra/src \
		-I$(top_srcdir)/../triutils/src \
		-I$(top_builddir)/../triutils/src

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libbelosutils_a_AR = $(ALTERNATE_AR)
else

libbelosutils_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# On alpha-dec machines, the extra .o's are needed for select packages only.
# ------------------------------------------------------------------------

libbelosutils_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(BELOSUTILS_LIB)

