
# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package TriKota subdirectory src built successfully."
	@echo ""

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

TRIKOTA_SRC_HPP = \
  $(srcdir)/TriKota_Version.hpp \
  $(srcdir)/TriKota_Driver.hpp \
  $(srcdir)/TriKota_DirectApplicInterface.hpp

TRIKOTA_SRC = \
  $(srcdir)/TriKota_Version.cpp \
  $(srcdir)/TriKota_Driver.cpp \
  $(srcdir)/TriKota_DirectApplicInterface.cpp


#Specify locations of include files within Trilinos

include $(top_builddir)/Makefile.export.TriKota

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(TRIKOTA_INCLUDES))
else
EXPORT_INCLUDES = $(TRIKOTA_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# --------------------------------------------------------------------
# TriKota library specifications
# --------------------------------------------------------------------

TRIKOTA_LIB = libtrikota.a

TRIKOTA_H = \
  $(srcdir)/TriKota_ConfigDefs.hpp \
  $(TRIKOTA_SRC_HPP) \
  $(TRIKOTA_SRC_H)

libtrikota_a_SOURCES = \
  $(TRIKOTA_SRC)

#We need the following to support separate build directories

# --------------------------------------------------------------------
# For using a special archiver
# --------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libtrikota_a_AR = $(ALTERNATE_AR)

else

libtrikota_a_AR = $(AR) $(ARFLAGS)

endif

# --------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# --------------------------------------------------------------------

# I have not been able to find any kind of "or" for Automake
# conditionals. JW

 libtrikota_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# --------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# --------------------------------------------------------------------

lib_LIBRARIES = $(TRIKOTA_LIB)

# --------------------------------------------------------------------
# List of all headers to install in $(includedir)
# --------------------------------------------------------------------

include_HEADERS = $(TRIKOTA_H) 

# --------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# --------------------------------------------------------------------

# SPECIAL NOTE: TriKotaConfig.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

execincludedir = $(exec_prefix)/include
nodist_execinclude_HEADERS = TriKota_config.h

# --------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# --------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 
