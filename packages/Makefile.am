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

if BUILD_AMESOS
AMESOS_SUBDIR = amesos
AMESOS_LIB_BUILD_CMND = cd amesos ; $(MAKE) all-libs ;
endif

if BUILD_ANASAZI
ANASAZI_SUBDIR = anasazi
ANASAZI_LIB_BUILD_CMND = cd anasazi ; $(MAKE) all-libs ;
endif

if BUILD_AZTECOO
AZTECOO_SUBDIR = aztecoo
AZTECOO_LIB_BUILD_CMND = cd aztecoo ; $(MAKE) all-libs ;
endif

if BUILD_BELOS
BELOS_SUBDIR = belos
BELOS_LIB_BUILD_CMND = cd belos ; $(MAKE) all-libs ;
endif

if BUILD_CAPO
CAPO_SUBDIR = capo
CAPO_LIB_BUILD_CMND = echo "capo does not support all-libs yet";
endif

if BUILD_CLAPS
CLAPS_SUBDIR = claps
CLAPS_LIB_BUILD_CMND = echo "claps does not support all-libs yet";
endif

if BUILD_DIDASKO
DIDASKO_SUBDIR = didasko
DIDASKO_LIB_BUILD_CMND = echo "disasko does not support all-libs yet";
endif

if BUILD_EPETRA
EPETRA_SUBDIR = epetra
EPETRA_LIB_BUILD_CMND = cd epetra ; $(MAKE) all-libs ;
endif

if BUILD_EPETRAEXT
EPETRAEXT_SUBDIR = epetraext
EPETRAEXT_LIB_BUILD_CMND = cd epetraext ; $(MAKE) all-libs ;
endif

if BUILD_IFPACK
IFPACK_SUBDIR = ifpack
IFPACK_LIB_BUILD_CMND = cd ifpack ; $(MAKE) all-libs ;
endif

if BUILD_KOKKOS
KOKKOS_SUBDIR = kokkos
KOKKOS_LIB_BUILD_CMND = echo "kokkos does not support all-libs yet";
endif

if BUILD_KOMPLEX
KOMPLEX_SUBDIR = komplex
KOMPLEX_LIB_BUILD_CMND = echo "komplex does not support all-libs yet";
endif

if BUILD_MEROS
MEROS_SUBDIR = meros
MEROS_LIB_BUILD_CMND = echo "meros does not support all-libs yet";
endif

if BUILD_ML
ML_SUBDIR = ml
ML_LIB_BUILD_CMND = echo "ml does not support all-libs yet";
endif

if BUILD_NEW_PACKAGE
NEW_PACKAGE_SUBDIR = new_package
NEW_PACKAGE_LIB_BUILD_CMND = echo "new_package does not support all-libs yet";
endif

if BUILD_NOX
NOX_SUBDIR = nox
NOX_LIB_BUILD_CMND = cd nox ; $(MAKE) all-libs ;
endif

if BUILD_PLIRIS
PLIRIS_SUBDIR = pliris
PLIRIS_LIB_BUILD_CMND = echo "pliris does not support all-libs yet";
endif

if BUILD_PYTRILINOS
PYTRILINOS_SUBDIR = PyTrilinos 
PYTRILINOS_LIB_BUILD_CMND = cd PyTrilinos ; $(MAKE) all-libs ;
endif

if BUILD_RYTHMOS
RYTHMOS_SUBDIR = rythmos
RYTHMOS_LIB_BUILD_CMND = cd rythmos ; $(MAKE) all-libs ;
endif

if BUILD_TEUCHOS
TEUCHOS_SUBDIR = teuchos
TEUCHOS_LIB_BUILD_CMND = cd teuchos ; $(MAKE) all-libs ;
endif

if BUILD_THYRA
THYRA_SUBDIR = thyra
THYRA_LIB_BUILD_CMND = cd thyra ; $(MAKE) all-libs ;
endif

if BUILD_TPETRA
TPETRA_SUBDIR = tpetra
TPETRA_LIB_BUILD_CMND = echo "tpetra does not support all-libs yet";
endif

if BUILD_TRIUTILS
TRIUTILS_SUBDIR = triutils
TRIUTILS_LIB_BUILD_CMND = cd triutils ; $(MAKE) all-libs ;
endif

if BUILD_TSF
TSF_SUBDIR = TSF
TSF_LIB_BUILD_CMND = echo "TSF does not support all-libs yet";
endif

if BUILD_TSFCORE
TSFCORE_SUBDIR = TSFCore
TSFCORE_LIB_BUILD_CMND = echo "TSFCore does not support all-libs yet";
endif

if BUILD_TSFCOREUTILS
TSFCOREUTILS_SUBDIR = TSFCoreUtils
TSFCOREUTILS_LIB_BUILD_CMND = echo "TSFCoreUtils does not support all-libs yet";
endif

if BUILD_TSFEXTENDED
TSFEXTENDED_SUBDIR = TSFExtended
TSFEXTENDED_LIB_BUILD_CMND = echo "TSFExtended does not support all-libs yet";
endif

# Must list the packages in the order that they are built!

SUBDIRS = \
	$(TEUCHOS_SUBDIR)      \
	$(PYTRILINOS_SUBDIR)   \
	$(THYRA_SUBDIR)        \
	$(KOKKOS_SUBDIR)       \
	$(EPETRA_SUBDIR)       \
	$(TRIUTILS_SUBDIR)     \
	$(TPETRA_SUBDIR)       \
	$(EPETRAEXT_SUBDIR)    \
	$(AZTECOO_SUBDIR)      \
	$(AMESOS_SUBDIR)       \
	$(IFPACK_SUBDIR)       \
	$(KOMPLEX_SUBDIR)      \
	$(CLAPS_SUBDIR)        \
	$(TSFCOREUTILS_SUBDIR) \
	$(TSFCORE_SUBDIR)      \
	$(BELOS_SUBDIR)        \
	$(ANASAZI_SUBDIR)      \
	$(PLIRIS_SUBDIR)       \
	$(ML_SUBDIR)           \
	$(TSF_SUBDIR)          \
	$(MEROS_SUBDIR)        \
	$(NOX_SUBDIR)          \
	$(CAPO_SUBDIR)         \
	$(RYTHMOS_SUBDIR)      \
	$(NEW_PACKAGE_SUBDIR)  \
	$(TSFEXTENDED_SUBDIR)  \
	$(DIDASKO_SUBDIR)

all-libs:
	$(TEUCHOS_LIB_BUILD_CMND)
	$(PYTRILINOS_LIB_BUILD_CMND)
	$(THYRA_LIB_BUILD_CMND)
	$(KOKKOS_LIB_BUILD_CMND)
	$(EPETRA_LIB_BUILD_CMND)
	$(TRIUTILS_LIB_BUILD_CMND)
	$(TPETRA_LIB_BUILD_CMND)
	$(EPETRAEXT_LIB_BUILD_CMND)
	$(AZTECOO_LIB_BUILD_CMND)
	$(AMESOS_LIB_BUILD_CMND)
	$(IFPACK_LIB_BUILD_CMND)
	$(KOMPLEX_LIB_BUILD_CMND)
	$(CLAPS_LIB_BUILD_CMND)
	$(TSFCOREUTILS_LIB_BUILD_CMND)
	$(TSFCORE_LIB_BUILD_CMND)
	$(BELOS_LIB_BUILD_CMND)
	$(ANASAZI_LIB_BUILD_CMND)
	$(PLIRIS_LIB_BUILD_CMND)
	$(ML_LIB_BUILD_CMND)
	$(TSF_LIB_BUILD_CMND)
	$(MEROS_LIB_BUILD_CMND)
	$(NOX_LIB_BUILD_CMND)
	$(CAPO_LIB_BUILD_CMND)
	$(RYTHMOS_LIB_BUILD_CMND)
	$(NEW_PACKAGE_LIB_BUILD_CMND)
	$(TSFEXTENDED_LIB_BUILD_CMND)
	$(DIDASKO_LIB_BUILD_CMND)
