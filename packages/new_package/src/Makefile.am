# @HEADER
# @HEADER


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

#
# There is, at present, no distinction between most of the categories, 
# other than those which start with EXTRA_.  All categories which do not
# start with EXTRA_ are included in all builds.  Categories which start
# with EXTRA_ are included in some builds (typically controlled by 
# --with-xyz or --without-xyz)
#
CORE_H = \
	$(srcdir)/Newp_Hello.h 

CORE = \
	$(srcdir)/Newp_Hello.cpp

UTIL_H = 

UTIL = 

ENUM_H = 

UNIPROC_H = 

UNIPROC = 

SERIAL_H = 

SERIAL = 

#
#  Files which are included in select builds are typically have
#  names that start with EXTRA
#
#  EXTRA_MPI is control by --with-mpi and --without-mpi which is
#  set up by TAC_ARC_CONFIG_MPI 
#
EXTRA_MPI_H =

EXTRA_MPI = 

if HAVE_MPI 
MPI_H = \
	$(EXTRA_MPI_H)
MPI = \
	$(EXTRA_MPI)

else

MPI_H =
MPI =

endif

EXTRA_NEWP_SWAHILI_H = \
	$(srcdir)/Newp_Jambo.h 

EXTRA_NEWP_SWAHILI = \
	$(srcdir)/Newp_Jambo.cpp

if HAVE_NEWP_SWAHILI

NEWP_SWAHILI_H = \
	$(EXTRA_NEWP_SWAHILI_H)
NEWP_SWAHILI = \
	$(EXTRA_NEWP_SWAHILI)

else

NEWP_SWAHILI_H =
NEWP_SWAHILI =

endif




USERBASE_H = 


LANG_WRAPS_H = 

LANG_WRAPS = 

DIST_H = 

DIST = 

COMP_H =

COMP = 


#Specify locations of include files within Trilinos, but outside New_Package
# #np# - Add packages that new_package is dependent on
#        Remove epetra if new_package does not need it
#  
AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src
# ------------------------------------------------------------------------
# New_Package library specifications
# ------------------------------------------------------------------------

NEW_PACKAGE_LIB = libnew_package.a


NEW_PACKAGE_H = \
	$(CORE_H) \
	$(UTIL_H) \
	$(ENUM_H) \
	$(UNIPROC_H) \
	$(SERIAL_H) \
	$(USERBASE_H) \
	$(DIST_H) \
	$(COMP_H) \
	$(MPI_H) \
	$(NEWP_SWAHILI_H)

libnew_package_a_SOURCES = \
	$(CORE) \
	$(UTIL) \
	$(UNIPROC) \
	$(SERIAL) \
	$(DIST) \
	$(COMP) \
	$(MPI) \
	$(NEWP_SWAHILI)

EXTRA_libnew_package_a_SOURCES = \
	$(EXTRA_MPI) \
	$(EXTRA_MPI_H) \
	$(EXTRA_NEWP_SWAHILI) \
	$(EXTRA_NEWP_SWAHILI_H)

#We need the following to support separate build directories
#AM_CPPFLAGS=-I$(builddir)/src
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libnew_package_a_AR = $(ALTERNATE_AR)
else

libnew_package_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libnew_package_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(NEW_PACKAGE_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(pkgincludedir)
# ------------------------------------------------------------------------

pkginclude_HEADERS = $(NEW_PACKAGE_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(pkgexecincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: New_Package_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

pkgexecincludedir = $(exec_prefix)/include/new_package
nodist_pkgexecinclude_HEADERS = New_Package_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






