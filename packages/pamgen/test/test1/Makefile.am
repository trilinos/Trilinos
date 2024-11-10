
include $(top_builddir)/Makefile.export.pamgen

# The following line helps the test harness recover from build errors.

all-local: 

EXEEXT = .exe

noinst_PROGRAMS = pamgen_test_driver 

pamgen_test_driver_SOURCES = $(srcdir)/../../example/pamgen_test_driver.c $(srcdir)/../../example/getopts.C

#np# Specify locations of include files within Trilinos, but outside np
#np# Add packages that np is dependent on; remove epetra if np does not need it

AM_CPPFLAGS = $(PAMGEN_INCLUDES)

#np# - List libraries that np is dependent on.  The DEPENDENCIES primary is
#np# optional.
pamgen_test_driver_DEPENDENCIES = \
	$(top_builddir)/src/libpamgen.a 

#np# List libraries that np is dependent on.
pamgen_test_driver_LDADD = $(PAMGEN_LIBS)

#
#  This is a crude but effective way to move files from the source tree to the build tree
#  We recommend renaming files because, although we recommend against the practice, some
#  users will configure and build in the same tree.
#
#  This requires at least one file to end in .exe (or whatever EXEEXT is set to).
#  The other files can have arbitrary names. (In this case Suppressions)
#
#np# Set "SUBDIRS= (list all subdirectories that are autotool'ed here)"
#np# List all subdirectories even if one or more are configured/built only
#np# conditionally.
#np# For packages whose examples are located in this directory, rather than
#np# in a subdirectory of this directory, use the Makefile.am located in
#np# example1/ rather than this file.
# SUBDIRS = example1
SUBDIRS = 
