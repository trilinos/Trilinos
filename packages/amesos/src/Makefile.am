# @HEADER
# @HEADER


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside Amesos. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

if HAVE_AMESOS_SLUD
SLUDIST_H = \
        $(srcdir)/SuperludistOO.h
SLUDIST = \
        $(srcdir)/SuperludistOO.cpp
else
SLUDIST_H =
SLUDIST =
endif

if BUILD_SLUSERIAL
SLUSERIAL_H = \
       $(srcdir)/Epetra_SLU.h
SLUSERIAL = \
       $(srcdir)/Epetra_SLU.cpp
else
SLUSERIAL_H =
SLUSERIAL =
endif

if BUILD_KUNDERT
KUNDERT_H = \
       $(srcdir)/Epetra_CrsKundertSparse.h
KUNDERT = \
       $(srcdir)/Epetra_CrsKundertSparse.cpp \
       $(srcdir)/ksparse/alloc.c \
       $(srcdir)/ksparse/smplinkrows.c \
       $(srcdir)/ksparse/spbuild.c \
       $(srcdir)/ksparse/spoutput.c \
       $(srcdir)/ksparse/spsolve.c \
       $(srcdir)/ksparse/driver.c \
       $(srcdir)/ksparse/spalloc.c \
       $(srcdir)/ksparse/spfactor.c \
       $(srcdir)/ksparse/spsmp.c \
       $(srcdir)/ksparse/sputils.c
else
KUNDERT_H =
KUNDERT =
endif

COMP_H = \
	$(srcdir)/CrsMatrixTranspose.h
COMP = \
	$(srcdir)/CrsMatrixTranspose.cpp

#
#  I turned off the addition of external/slu to the include path 
#  so that I could get SLU dist to work.  Ken Feb 11, 2003
#
#Specify locations of include files within Trilinos, but outside Amesos
# AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src \
#               -I$(top_srcdir)/../../external/slu -I$(top_srcdir)/src/ksparse

if BUILD_SLUSERIAL
SLUSERIAL_CPPFLAGS = -I$(top_srcdir)/../../external/slu
else
SLUSERIAL_CPPFLAGS = 
endif

if BUILD_KUNDERT
KUNDERT_CPPFLAGS = -I$(top_srcdir)/src/ksparse
else
KUNDERT_CPPFLAGS = 
endif

AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src \
              $(SLUSERIAL_CPPFLAGS) $(KUNDERT_CPPFLAGS)

# ------------------------------------------------------------------------
# Amesos library specifications
# ------------------------------------------------------------------------

AMESOS_LIB = libamesos.a


AMESOS_H = \
        $(SLUDIST_H) \
        $(SLUSERIAL_H) \
        $(KUNDERT_H) \
	$(COMP_H)

libamesos_a_SOURCES = \
        $(SLUDIST) \
        $(SLUSERIAL) \
        $(KUNDERT) \
	$(COMP)

EXTRA_libamesos_a_SOURCES = \
	$(EXTRA_MPI) \
	$(EXTRA_MPI_H) \
	$(EXTRA_MPISMP) \
	$(EXTRA_MPISMP_H)

#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libamesos_a_AR = $(ALTERNATE_AR)
else

libamesos_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libamesos_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(AMESOS_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(pkgincludedir)
# ------------------------------------------------------------------------

pkginclude_HEADERS = $(AMESOS_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(pkgexecincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: Amesos_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

pkgexecincludedir = $(exec_prefix)/include/amesos
nodist_pkgexecinclude_HEADERS = Amesos_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






