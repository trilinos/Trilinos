# @HEADER
# @HEADER

SOURCE_H = \
	$(srcdir)/Thyra_MLPreconditionerFactory.hpp 

SOURCE = \
	$(srcdir)/Thyra_MLPreconditionerFactory.cpp 

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
# Ml Thyra library specifications
# --------------------------------------------------------------------

ML_THYRA_LIB = libstratimikosml.a

ML_THYRA_H = \
	$(SOURCE_H)

libstratimikosml_a_SOURCES = \
	$(SOURCE)

#We need the following to support separate build directories

# --------------------------------------------------------------------
# For using a special archiver
# --------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libstratimikosml_a_AR = $(ALTERNATE_AR)

else

libstratimikosml_a_AR = $(AR) cru

endif

# --------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# --------------------------------------------------------------------

# I have not been able to find any kind of "or" for Automake
# conditionals. JW

 libstratimikosml_a_LIBADD = $(TRILINOS_TEMPLATE_OBJS) $(XTRALDADD)

# --------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# --------------------------------------------------------------------

lib_LIBRARIES = $(ML_THYRA_LIB)

# --------------------------------------------------------------------
# List of all headers to install in $(includedir)
# --------------------------------------------------------------------

include_HEADERS = $(ML_THYRA_H) 

# --------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# --------------------------------------------------------------------

execincludedir = $(exec_prefix)/include

# --------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# --------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 
