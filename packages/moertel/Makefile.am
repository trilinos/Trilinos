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
ACLOCAL_AMFLAGS = -I config

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
config/replace-install-prefix.pl 	  config/string-replace.pl \
config/tac_arg_enable_export-makefiles.m4 config/token-replace.pl \
config/ac_cxx_std_sprintf.m4 \
config/tac_arg_with_ar.m4 \
config/tac_arg_with_flags.m4 \
config/tac_arg_with_incdirs.m4 \
config/config.guess                       config/tac_arg_with_libdirs.m4 \
config/config.sub \
config/strip_dup_incl_paths.pl            config/strip_dup_libs.pl \
config/tac_arg_with_perl.m4 \
bootstrap                                 doc/build_docs \
doc/Doxyfile                              doc/index.doc
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
	@if [ $(srcdir) != '.' ]; then cp -f $(srcdir)/example/TwoSquares.grid $(top_builddir)/example/TwoSquares.grid; cp -f $(srcdir)/example/TwoSquares.m $(top_builddir)/example/TwoSquares.m; cp -f $(srcdir)/example/Hex_3D.grid $(top_builddir)/example/Hex_3D.grid; fi
	@echo "        +----------------------------------------------------+"
	@echo "        |    Trilinos package Moertel built successfully     |"
	@echo "        +----------------------------------------------------+"

## #######################################################################
## Export Makefile Installation
## #######################################################################
if USING_EXPORT_MAKEFILES

install-exec-hook:
	mkdir -p $(DESTDIR)$(includedir)
	cp $(top_builddir)/Makefile.export.moertel $(DESTDIR)$(includedir)/.
	$(PERL_EXE) $(top_srcdir)/config/replace-install-prefix.pl \
		--exec-prefix=$(exec_prefix) \
		--my-export-makefile=Makefile.export.moertel \
		--my-abs-top-srcdir=@abs_top_srcdir@ \
		--my-abs-incl-dirs=@abs_top_builddir@/src:@abs_top_srcdir@/src \
		--my-abs-lib-dirs=@abs_top_builddir@/src \
		--dep-package-abs-builddirs=@abs_top_builddir@/../epetra:@abs_top_builddir@/../teuchos:@abs_top_builddir@/../epetraext:@abs_top_builddir@/../ml:@abs_top_builddir@/../amesos
	$(PERL_EXE) $(top_srcdir)/config/generate-makeoptions.pl $(top_builddir)/src/Makefile \
		MOERTEL > $(DESTDIR)$(includedir)/Makefile.export.moertel.macros

uninstall-hook:
	rm -f $(includedir)/Makefile.export.moertel
	rm -f $(includedir)/Makefile.export.moertel.macros

else

install-exec-hook:

uninstall-hook:

endif

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################

#We now build tests and examples through separate make targets, rather than
#during "make".  We still need to conditionally include the test and example
#in SUBDIRS, even though BUILD_TESTS and BUILD_EXAMPLES will never be
#defined, so that the tests and examples are included in the distribution
#tarball.
#(Moertel does not yet have tests.)

#SUB_TEST
#TEST_SUBDIR=test
#endif

if SUB_EXAMPLE
EXAMPLE_SUBDIR=example
endif

tests:
	@echo "moertel does not have any tests yet"

if BUILD_EXAMPLES
examples:
	@echo ""
	@echo "Now building moertel examples."
	@echo ""
	cd $(top_builddir)/example && $(MAKE)
	@echo ""
	@echo "Finished building moertel examples."
	@echo ""

install-examples:
	cd $(top_builddir)/example && $(MAKE) install
else
examples:
	@echo "moertel examples were disabled at configure time"

install-examples:
	@echo "moertel examples were disabled at configure time"
endif

clean-tests:
	@echo "moertel does not have any tests yet"

clean-examples:
	cd $(top_builddir)/example && $(MAKE) clean

everything:
	$(MAKE) && $(MAKE) examples && $(MAKE) tests

clean-everything:
	$(MAKE) clean-examples && $(MAKE) clean-tests && $(MAKE) clean

install-everything:
	$(MAKE) install && $(MAKE) install-examples

SUBDIRS = src $(EXAMPLE_SUBDIR)

TRILINOS_HOME_DIR=@abs_top_srcdir@/../..
TRILINOS_BUILD_DIR=@abs_top_builddir@/../..

