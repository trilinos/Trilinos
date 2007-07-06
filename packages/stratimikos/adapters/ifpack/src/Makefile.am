
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

SOURCE_H = \
	$(srcdir)/Thyra_IfpackPreconditionerFactory.hpp

SOURCE = \
	$(srcdir)/Thyra_IfpackPreconditionerFactory.cpp

# Specify locations of include files within Trilinos

include $(top_builddir)/Makefile.export.stratimikos

_EXPORT_INCLUDES=$(STRATIMIKOS_INCLUDES)

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(_EXPORT_INCLUDES))
else
EXPORT_INCLUDES = $(_EXPORT_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# --------------------------------------------------------------------
# Ifpack Thyra library specifications
# --------------------------------------------------------------------

IFPACK_THYRA_LIB = libstratimikosifpack.a

IFPACK_THYRA_H = \
	$(top_srcdir)/src/Stratimikos_ConfigDefs.hpp \
	$(SOURCE_H)

libstratimikosifpack_a_SOURCES = \
	$(SOURCE)

#We need the following to support separate build directories

# --------------------------------------------------------------------
# For using a special archiver
# --------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libstratimikosifpack_a_AR = $(ALTERNATE_AR)

else

libstratimikosifpack_a_AR = $(AR) cru

endif

# --------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# --------------------------------------------------------------------

# I have not been able to find any kind of "or" for Automake
# conditionals. JW

 libstratimikosifpack_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# --------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# --------------------------------------------------------------------

lib_LIBRARIES = $(IFPACK_THYRA_LIB)

# --------------------------------------------------------------------
# List of all headers to install in $(includedir)
# --------------------------------------------------------------------

include_HEADERS = $(IFPACK_THYRA_H) 

# --------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# --------------------------------------------------------------------

execincludedir = $(exec_prefix)/include

# --------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# --------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 
