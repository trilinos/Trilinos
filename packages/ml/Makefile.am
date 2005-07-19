# Typically the only line that needs to be edited in this file is the last 
# line.  Set 'SUBDIRS= (list all subdirectories that are autotool'ed here)'
# List all subdirectories even if one or more are configured/built only
# conditionally.

# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001,2002,2003) Sandia Corporation
#

# Copyright (2001,2002,2003) Sandia Corportation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
#
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in ths data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from

# others acting on its behalf a paid-up, nonexclusive, irrevocable
# worldwide license in this data to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
#
# NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
# OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
# ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
# RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
# INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
# THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
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
EXTRA_DIST = bootstrap \
	config/ac_cxx_have_std.m4 config/ac_cxx_have_stl.m4 \
	config/ac_cxx_namespaces.m4 config/ac_cxx_mutable.m4 \
	config/depcomp ChangeLog-ML README-ML COPYRIGHT

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
	@echo "Trilinos package ml built successfully"

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################

if BUILD_TESTS
TEST_SUBDIR=test
endif

if ML_BUILD_EXAMPLES
EXAMPLE_SUBDIR=examples
endif

if PYTHON_USE
PYTHON_SUBDIR=python
endif

SUBDIRS = src $(EXAMPLE_SUBDIR) $(TEST_SUBDIR) doc $(PYTHON_SUBDIR)


## #######################################################################
## Target for circular dependencies
## #######################################################################
all-libs:
	cd src ; make all-libs

## #######################################################################
## Export Makefile Installation
## #######################################################################
if USING_EXPORT_MAKEFILES

install-exec-hook:
	cp $(top_builddir)/Makefile.export.ml $(exec_prefix)/include/.
	$(PERL_EXE) $(top_srcdir)/config/replace-install-prefix.pl \
		--exec-prefix=$(exec_prefix) \
		--my-export-makefile=Makefile.export.ml \
		--my-abs-top-srcdir=@abs_top_srcdir@ \
		--my-abs-incl-dirs=@abs_top_builddir@/src:@abs_top_srcdir@/src:@abs_top_srcdir@/Comm:abs_top_srcdir@/Operator:abs_top_srcdir@/Smoother:@abs_top_srcdir@/Coarsen:@abs_top_srcdir@/Krylov:@abs_top_srcdir@/Main:@abs_top_srcdir@/FEGrid:@abs_top_srcdir@/Utils:@abs_top_srcdir@/MLAPI:@abs_top_srcdir@/Include \
		--my-abs-lib-dirs=@abs_top_builddir@/src \
		--dep-package-abs-builddirs=@abs_top_builddir@/../epetraext:@abs_top_builddir@/../epetra:@abs_top_builddir@/../triutils:@abs_top_builddir@/../amesos:@abs_top_builddir@/../ifpack:@abs_top_builddir@/../teuchos:@abs_top_builddir@/../anasazi:@abs_top_builddir@/../aztecoo
	$(PERL_EXE) $(top_srcdir)/config/generate-makeoptions.pl $(top_builddir)/src/Makefile ml > $(exec_prefix)/include/Makefile.export.ml.macros

uninstall-hook:
	rm -f $(exec_prefix)/include/Makefile.export.ml
	rm -f $(exec_prefix)/include/Makefile.export.ml.macros

else 

install-exec-hook:

uninstall-hook:

endif
