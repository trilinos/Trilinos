# @HEADER
# ************************************************************************
# 
#                          Moertel FE Package
#                 Copyright (2006) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael Gee (mwgee@sandia.gov) 
# 
# ************************************************************************
# @HEADER

## #######################################################################
## Options to automake (rarely used - don't worry about it)
## #######################################################################
AUTOMAKE_OPTIONS = foreign

## #######################################################################
## Aclocal command (rarely used - don't worry about it)
## #######################################################################
ACLOCAL	= aclocal -I config

#
#       I believe that by switching to AUX_DIR(../../config) one 
#       could get rid of these.
#
## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
#np# For a typical package, there is no reason to distribute these files
#np# because users should not have to bootstrap.  We distribute them with
#np# new package so that the files can be used in creating the 
#np# configure script for other packages. 
EXTRA_DIST = \
config/ac_cxx_bool.m4                     config/depcomp \
config/install-sh \
config/ac_cxx_have_stl.m4                 config/missing \
config/ac_cxx_mutable.m4                  config/mkinstalldirs \
config/ac_cxx_namespaces.m4               config/tac_arg_check_mpi.m4 \
config/tac_arg_config_mpi.m4 		  config/tac_arg_enable_feature.m4 \
config/generate-makeoptions.pl    config/tac_arg_enable_feature_sub_check.m4 \
config/replace-install-prefix.pl \
config/tac_arg_enable_export-makefiles.m4 \
config/ac_cxx_std_sprintf.m4 \
config/tac_arg_with_ar.m4 \
config/tac_arg_with_flags.m4 \
config/tac_arg_with_incdirs.m4 \
config/config.guess                       config/tac_arg_with_libdirs.m4 \
config/config.sub \
config/strip_dup_incl_paths.pl            config/strip_dup_libs_.pl \
config/tac_arg_with_perl.pl \
bootstrap 
#
#       I believe that by switching to AUX_DIR(../../config) one 
#       could get rid of these.
#
## #######################################################################
## Tools in the auxillary directory 
## #######################################################################
AUX_DIST = config/install-sh config/missing config/mkinstalldirs 
#
#  Again, I hope that AUX_DIR(../../config) eliminates these
#  config/install-sh config/missing config/mkinstalldirs 

## #######################################################################
## Files to be deleted by 'make maintainer-clean'
## #######################################################################
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 autom4te.cache/* \
	configure config.status config.log \
	src/common/config-h.in src/common/stamp-h.in \
	$(AUX_DIST) 

#The following line helps the test harness recover from build errors.
                                                                                
all-local:
	@cp -f $(srcdir)/example/TwoSquares.grid $(top_builddir)/example/TwoSquares.grid
	@cp -f $(srcdir)/example/TwoSquares.m $(top_builddir)/example/TwoSquares.m
	@cp -f $(srcdir)/example/Hex_3D.grid $(top_builddir)/example/Hex_3D.grid
	@echo "        +----------------------------------------------------+"
	@echo "        |    Trilinos package Moertel built successfully     |"
	@echo "        +----------------------------------------------------+"

## #######################################################################
## Target for circular dependencies
## #######################################################################
all-libs:
	cd src ; make all-libs

libs:
	cd src ; make

## #######################################################################
## Export Makefile Installation
## #######################################################################
if USING_EXPORT_MAKEFILES

install-exec-hook:
	cp $(top_builddir)/Makefile.export.moertel $(exec_prefix)/include/.
	$(PERL_EXE) $(top_srcdir)/config/replace-install-prefix.pl \
		--exec-prefix=$(exec_prefix) \
		--my-export-makefile=Makefile.export.moertel \
		--my-abs-top-srcdir=@abs_top_srcdir@ \
		--my-abs-incl-dirs=@abs_top_builddir@/src:@abs_top_srcdir@/src \
		--my-abs-lib-dirs=@abs_top_builddir@/src \
		--dep-package-abs-builddirs=@abs_top_builddir@/../epetra:@abs_top_builddir@/../teuchos:@abs_top_builddir@/../epetraext:@abs_top_builddir@/../ml:@abs_top_builddir@/../amesos
	$(PERL_EXE) $(top_srcdir)/config/generate-makeoptions.pl $(top_builddir)/src/Makefile \
		MOERTEL > $(exec_prefix)/include/Makefile.export.moertel.macros

uninstall-hook:
	rm -f $(exec_prefix)/include/Makefile.export.moertel
	rm -f $(exec_prefix)/include/Makefile.export.moertel.macros

else

install-exec-hook:

uninstall-hook:

endif

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################

#if BUILD_TESTS
#TEST_SUBDIR=test
#endif

#if BUILD_EXAMPLES
#EXAMPLE_SUBDIR=example
#endif

#if PYTHON_USE
#PYTHON_SUBDIR=python
#endif

#SUBDIRS = src $(TEST_SUBDIR) $(EXAMPLE_SUBDIR) $(PYTHON_SUBDIR)
SUBDIRS = src

TRILINOS_HOME_DIR=@abs_top_srcdir@/../..
TRILINOS_BUILD_DIR=@abs_top_builddir@/../..

