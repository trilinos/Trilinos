# @HEADER
# ************************************************************************
#
#            Trilinos: An Object-Oriented Solver Framework
#                 Copyright (2001) Sandia Corporation
#

# Copyright (2001) Sandia Corportation. Under the terms of Contract
# DE-AC04-94AL85000, there is a non-exclusive license for use of this
# work by or on behalf of the U.S. Government.  Export of this program
# may require a license from the United States Government.
#
# NOTICE:  The United States Government is granted for itself and others
# acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
# license in ths data to reproduce, prepare derivative works, and
# perform publicly and display publicly.  Beginning five (5) years from
# July 25, 2001, the United States Government is granted for itself and
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

## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
#EXTRA_DIST = bootstrap \
#	config/ac_cxx_have_std.m4 config/ac_cxx_have_stl.m4 \
#	config/ac_cxx_namespaces.m4 config/ac_cxx_mutable.m4 

## #######################################################################
## Tools in the auxillary directory 
## #######################################################################
#AUX_DIST = config/install-sh config/missing config/mkinstalldirs

## #######################################################################
## Files to be deleted by 'make maintainer-clean'
## #######################################################################
#MAINTAINERCLEANFILES = Makefile.in aclocal.m4 autom4te.cache/* \
#	configure config.status config.log \
#	src/common/config-h.in src/common/stamp-h.in \
#	$(AUX_DIST) 

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################
#SUBDIRS = aztecoo epetra ifpack komplex petra tpetra triutils tsf TSF epetraext

if BUILD_KOKKOS
KOKKOS_SUBDIR = kokkos
endif

if BUILD_TSFCORE
TSFCORE_SUBDIR = TSFCore
endif

if BUILD_TSFEXTENDED
TSFEXTENDED_SUBDIR = TSFExtended
endif

if BUILD_TSFCOREUTILS
TSFCOREUTILS_SUBDIR = TSFCoreUtils
endif

if BUILD_EPETRAEXT
EPETRAEXT_SUBDIR = epetraext
endif

if BUILD_AMESOS
AMESOS_SUBDIR = amesos
endif

if BUILD_BELOS
BELOS_SUBDIR = belos
endif

if BUILD_ANASAZI
ANASAZI_SUBDIR = anasazi
endif

## #np# - Copy the next three lines, and change the name
if BUILD_NEW_PACKAGE
NEW_PACKAGE_SUBDIR = new_package
endif

if BUILD_TPETRA
TPETRA_SUBDIR = tpetra
endif

if BUILD_TEUCHOS
TEUCHOS_SUBDIR = teuchos
endif

if BUILD_NOX
NOX_SUBDIR = nox
endif

if BUILD_KOMPLEX
KOMPLEX_SUBDIR = komplex
endif

if BUILD_AZTECOO
AZTECOO_SUBDIR = aztecoo
endif

if BUILD_IFPACK
IFPACK_SUBDIR = ifpack
endif

if BUILD_TRIUTILS
TRIUTILS_SUBDIR = triutils
endif

if BUILD_EPETRA
EPETRA_SUBDIR = epetra
endif

if BUILD_TSF
TSF_SUBDIR = TSF
endif

if BUILD_MEROS
MEROS_SUBDIR = meros
endif

if BUILD_ML
ML_SUBDIR = ml
endif

#Epetra must be listed first due to aztecoo's dependence on Epetra
## #np# - Add your subdir, i.e. $(NEW_PACKAGE_SUBDIR) which is defined 
## above, to this list
SUBDIRS = $(TEUCHOS_SUBDIR) $(KOKKOS_SUBDIR) $(EPETRA_SUBDIR) \
	  $(TRIUTILS_SUBDIR) $(TPETRA_SUBDIR) $(IFPACK_SUBDIR) \
	  $(AZTECOO_SUBDIR) $(BELOS_SUBDIR) $(ANASAZI_SUBDIR) \
	  $(KOMPLEX_SUBDIR) \
	  $(TSFCOREUTILS_SUBDIR) $(TSFCORE_SUBDIR) $(TSF_SUBDIR) \
		$(TSFEXTENDED_SUBDIR) \
	  $(AMESOS_SUBDIR) $(ML_SUBDIR) $(EPETRAEXT_SUBDIR) \
	  $(MEROS_SUBDIR) $(NEW_PACKAGE_SUBDIR) \
          $(NOX_SUBDIR)


