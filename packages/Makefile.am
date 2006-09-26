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
AMESOS_LIBS_BUILD_CMND = cd amesos ; $(MAKE) libs ;
AMESOS_ALL_LIBS_BUILD_CMND = cd amesos ; $(MAKE) all-libs ;
endif

if BUILD_ANASAZI
ANASAZI_SUBDIR = anasazi
ANASAZI_LIBS_BUILD_CMND = cd anasazi ; $(MAKE) libs ;
ANASAZI_ALL_LIBS_BUILD_CMND = cd anasazi ; $(MAKE) all-libs ;
endif

if BUILD_AZTECOO
AZTECOO_SUBDIR = aztecoo
AZTECOO_LIBS_BUILD_CMND = cd aztecoo ; $(MAKE) libs ;
AZTECOO_ALL_LIBS_BUILD_CMND = cd aztecoo ; $(MAKE) all-libs ;
endif

if BUILD_BELOS
BELOS_SUBDIR = belos
BELOS_LIBS_BUILD_CMND = cd belos ; $(MAKE) libs ;
BELOS_ALL_LIBS_BUILD_CMND = cd belos ; $(MAKE) all-libs ;
endif

if BUILD_CAPO
CAPO_SUBDIR = capo
CAPO_LIBS_BUILD_CMND = echo "capo does not support libs target yet";
CAPO_ALL_LIBS_BUILD_CMND = echo "capo does not support all-libs target yet";
endif

if BUILD_CLAPS
CLAPS_SUBDIR = claps
CLAPS_LIBS_BUILD_CMND = echo "claps does not support libs target yet";
CLAPS_ALL_LIBS_BUILD_CMND = echo "claps does not support all-libs target yet";
endif

if BUILD_DIDASKO
DIDASKO_SUBDIR = didasko
DIDASKO_LIBS_BUILD_CMND = echo "disasko does not support all-libs target yet";
DIDASKO_ALL_LIBS_BUILD_CMND = echo "disasko does not support all-libs target yet";
endif

if BUILD_EPETRA
EPETRA_SUBDIR = epetra
EPETRA_LIBS_BUILD_CMND = cd epetra ; $(MAKE) libs ;
EPETRA_ALL_LIBS_BUILD_CMND = cd epetra ; $(MAKE) all-libs ;
endif

if BUILD_EPETRAEXT
EPETRAEXT_SUBDIR = epetraext
EPETRAEXT_LIBS_BUILD_CMND = cd epetraext ; $(MAKE) libs ;
EPETRAEXT_ALL_LIBS_BUILD_CMND = cd epetraext ; $(MAKE) all-libs ;
endif

if BUILD_EXTERNAL
EXTERNAL_SUBDIR = external
EXTERNAL_LIBS_BUILD_CMND = echo "external does not support libs target yet";
EXTERNAL_ALL_LIBS_BUILD_CMND = echo "external does not support all-libs target yet";
endif

if BUILD_IFPACK
IFPACK_SUBDIR = ifpack
IFPACK_LIBS_BUILD_CMND = cd ifpack ; $(MAKE) libs ;
IFPACK_ALL_LIBS_BUILD_CMND = cd ifpack ; $(MAKE) all-libs ;
endif

if BUILD_GALERI
GALERI_SUBDIR = galeri 
GALERI_LIBS_BUILD_CMND = cd galeri ; $(MAKE) libs ;
GALERI_ALL_LIBS_BUILD_CMND = cd galeri ; $(MAKE) all-libs ;
endif

if BUILD_KOKKOS
KOKKOS_SUBDIR = kokkos
KOKKOS_LIBS_BUILD_CMND = cd kokkos ; $(MAKE) libs ;
KOKKOS_ALL_LIBS_BUILD_CMND = cd kokkos ; $(MAKE) all-libs ;
endif

if BUILD_KOMPLEX
KOMPLEX_SUBDIR = komplex
KOMPLEX_LIBS_BUILD_CMND = echo "komplex does not support libs target yet";
KOMPLEX_ALL_LIBS_BUILD_CMND = echo "komplex does not support all-libs target yet";
endif

if BUILD_MEROS
MEROS_SUBDIR = meros
MEROS_LIBS_BUILD_CMND = echo "meros does not support libs target yet";
MEROS_ALL_LIBS_BUILD_CMND = echo "meros does not support all-libs target yet";
endif

if BUILD_ML
ML_SUBDIR = ml
ML_LIBS_BUILD_CMND = cd ml ; $(MAKE) libs ;
ML_ALL_LIBS_BUILD_CMND = cd ml ; $(MAKE) all-libs ;
endif

if BUILD_MOERTEL
MOERTEL_SUBDIR = moertel 
MOERTEL_LIBS_BUILD_CMND = cd moertel ; $(MAKE) libs ;
MOERTEL_ALL_LIBS_BUILD_CMND = cd moertel ; $(MAKE) all-libs ;
endif

if BUILD_MOOCHO
MOOCHO_SUBDIR = moocho
MOOCHO_LIBS_BUILD_CMND = cd moocho ; $(MAKE) libs ;
MOOCHO_ALL_LIBS_BUILD_CMND = cd moocho ; $(MAKE) all-libs ;
endif

if BUILD_NEW_PACKAGE
NEW_PACKAGE_SUBDIR = new_package
NEW_PACKAGE_LIBS_BUILD_CMND = cd new_package ; $(MAKE) libs ;
NEW_PACKAGE_ALL_LIBS_BUILD_CMND = cd new_package ; $(MAKE) all-libs ;
endif

if BUILD_ISORROPIA
ISORROPIA_SUBDIR = isorropia
ISORROPIA_ALL_LIBS_BUILD_CMND = cd isorropia ; $(MAKE) all-libs ;
endif

if BUILD_NOX
NOX_SUBDIR = nox
NOX_LIBS_BUILD_CMND = cd nox ; $(MAKE) libs ;
NOX_ALL_LIBS_BUILD_CMND = cd nox ; $(MAKE) all-libs ;
endif

if BUILD_PLIRIS
PLIRIS_SUBDIR = pliris
PLIRIS_LIBS_BUILD_CMND = echo "pliris does not support libs target yet";
PLIRIS_ALL_LIBS_BUILD_CMND = echo "pliris does not support all-libs target yet";
endif

if BUILD_PYTRILINOS
PYTRILINOS_SUBDIR = PyTrilinos 
PYTRILINOS_LIBS_BUILD_CMND = cd PyTrilinos ; $(MAKE) libs ;
PYTRILINOS_ALL_LIBS_BUILD_CMND = cd PyTrilinos ; $(MAKE) all-libs ;
endif

if BUILD_RBGEN
RBGEN_SUBDIR = rbgen
RBGEN_LIBS_BUILD_CMND = cd rbgen ; $(MAKE) libs ;
RBGEN_ALL_LIBS_BUILD_CMND = cd rbgen ; $(MAKE) all-libs ;
endif

if BUILD_RYTHMOS
RYTHMOS_SUBDIR = rythmos
RYTHMOS_LIBS_BUILD_CMND = cd rythmos ; $(MAKE) libs ;
RYTHMOS_ALL_LIBS_BUILD_CMND = cd rythmos ; $(MAKE) all-libs ;
endif

if BUILD_SACADO
SACADO_SUBDIR = sacado
SACADO_LIBS_BUILD_CMND = cd sacado ; $(MAKE) libs ;
SACADO_ALL_LIBS_BUILD_CMND = cd sacado ; $(MAKE) all-libs ;
endif

if BUILD_TEUCHOS
TEUCHOS_SUBDIR = teuchos
TEUCHOS_LIBS_BUILD_CMND = cd teuchos ; $(MAKE) libs ;
TEUCHOS_ALL_LIBS_BUILD_CMND = cd teuchos ; $(MAKE) all-libs ;
endif

if BUILD_RTOP
RTOP_SUBDIR = rtop
RTOP_LIBS_BUILD_CMND = cd rtop ; $(MAKE) libs ;
RTOP_ALL_LIBS_BUILD_CMND = cd rtop ; $(MAKE) all-libs ;
endif

if BUILD_THYRA
THYRA_SUBDIR = thyra
THYRA_LIBS_BUILD_CMND = cd thyra ; $(MAKE) libs ;
THYRA_ALL_LIBS_BUILD_CMND = cd thyra ; $(MAKE) all-libs ;
endif

if BUILD_STRATIMIKOS
STRATIMIKOS_SUBDIR = stratimikos
STRATIMIKOS_LIBS_BUILD_CMND = cd stratimikos ; $(MAKE) libs ;
STRATIMIKOS_ALL_LIBS_BUILD_CMND = cd stratimikos ; $(MAKE) all-libs ;
endif

if BUILD_TPETRA
TPETRA_SUBDIR = tpetra
TPETRA_LIBS_BUILD_CMND = cd tpetra ; $(MAKE) libs ;
TPETRA_ALL_LIBS_BUILD_CMND = cd tpetra ; $(MAKE) all-libs ;
endif

if BUILD_TRIUTILS
TRIUTILS_SUBDIR = triutils
TRIUTILS_LIBS_BUILD_CMND = cd triutils ; $(MAKE) libs ;
TRIUTILS_ALL_LIBS_BUILD_CMND = cd triutils ; $(MAKE) all-libs ;
endif

if BUILD_WEBTRILINOS
WEBTRILINOS_SUBDIR = WebTrilinos
WEBTRILINOS_LIBS_BUILD_CMND = echo "WebTrilinos does not support libs target";
WEBTRILINOS_ALL_LIBS_BUILD_CMND = echo "WebTrilinos does not support all-libs target";
endif

# Must list the packages in the order that they are built!

SUBDIRS = \
	$(TEUCHOS_SUBDIR)      \
	$(RTOP_SUBDIR)         \
	$(THYRA_SUBDIR)        \
	$(KOKKOS_SUBDIR)       \
	$(EPETRA_SUBDIR)       \
	$(TRIUTILS_SUBDIR)     \
	$(TPETRA_SUBDIR)       \
	$(EPETRAEXT_SUBDIR)    \
	$(ISORROPIA_SUBDIR)    \
	$(AZTECOO_SUBDIR)      \
	$(GALERI_SUBDIR)       \
	$(AMESOS_SUBDIR)       \
	$(IFPACK_SUBDIR)       \
	$(KOMPLEX_SUBDIR)      \
	$(CLAPS_SUBDIR)        \
	$(BELOS_SUBDIR)        \
	$(ANASAZI_SUBDIR)      \
        $(RBGEN_SUBDIR)        \
	$(PLIRIS_SUBDIR)       \
	$(ML_SUBDIR)           \
	$(MOERTEL_SUBDIR)      \
	$(STRATIMIKOS_SUBDIR)  \
	$(MEROS_SUBDIR)        \
	$(NOX_SUBDIR)          \
	$(CAPO_SUBDIR)         \
	$(RYTHMOS_SUBDIR)      \
	$(NEW_PACKAGE_SUBDIR)  \
	$(MOOCHO_SUBDIR)       \
	$(SACADO_SUBDIR)       \
	$(PYTRILINOS_SUBDIR)   \
	$(DIDASKO_SUBDIR) \
	$(WEBTRILINOS_SUBDIR) \
	$(EXTERNAL_SUBDIR)

libs:
	$(TEUCHOS_LIBS_BUILD_CMND)
	$(RTOP_LIBS_BUILD_CMND)
	$(THYRA_LIBS_BUILD_CMND)
	$(KOKKOS_LIBS_BUILD_CMND)
	$(EPETRA_LIBS_BUILD_CMND)
	$(TRIUTILS_LIBS_BUILD_CMND)
	$(TPETRA_LIBS_BUILD_CMND)
	$(EPETRAEXT_LIBS_BUILD_CMND)
	$(ISORROPIA_LIBS_BUILD_CMND)
	$(AZTECOO_LIBS_BUILD_CMND)
	$(AMESOS_LIBS_BUILD_CMND)
	$(IFPACK_LIBS_BUILD_CMND)
	$(GALERI_LIBS_BUILD_CMND)
	$(KOMPLEX_LIBS_BUILD_CMND)
	$(CLAPS_LIBS_BUILD_CMND)
	$(BELOS_LIBS_BUILD_CMND)
	$(ANASAZI_LIBS_BUILD_CMND)
	$(RBGEN_LIBS_BUILD_CMND)
	$(PLIRIS_LIBS_BUILD_CMND)
	$(ML_LIBS_BUILD_CMND)
	$(MOERTEL_LIBS_BUILD_CMND)
	$(STRATIMIKOS_LIBS_BUILD_CMND)
	$(MEROS_LIBS_BUILD_CMND)
	$(NOX_LIBS_BUILD_CMND)
	$(CAPO_LIBS_BUILD_CMND)
	$(RYTHMOS_LIBS_BUILD_CMND)
	$(NEW_PACKAGE_LIBS_BUILD_CMND)
	$(MOOCHO_LIBS_BUILD_CMND)
	$(SACADO_LIBS_BUILD_CMND)
	$(PYTRILINOS_LIBS_BUILD_CMND)
	$(DIDASKO_LIBS_BUILD_CMND)
	$(WEBTRILINOS_LIBS_BUILD_CMND)
	$(EXTERNAL_LIBS_BUILD_CMND)

all-libs:
	$(TEUCHOS_ALL_LIBS_BUILD_CMND)
	$(RTOP_ALL_LIBS_BUILD_CMND)
	$(THYRA_ALL_LIBS_BUILD_CMND)
	$(KOKKOS_ALL_LIBS_BUILD_CMND)
	$(EPETRA_ALL_LIBS_BUILD_CMND)
	$(TRIUTILS_ALL_LIBS_BUILD_CMND)
	$(TPETRA_ALL_LIBS_BUILD_CMND)
	$(EPETRAEXT_ALL_LIBS_BUILD_CMND)
	$(ISORROPIA_ALL_LIBS_BUILD_CMND)
	$(AZTECOO_ALL_LIBS_BUILD_CMND)
	$(AMESOS_ALL_LIBS_BUILD_CMND)
	$(IFPACK_ALL_LIBS_BUILD_CMND)
	$(GALERI_ALL_LIBS_BUILD_CMND)
	$(KOMPLEX_ALL_LIBS_BUILD_CMND)
	$(CLAPS_ALL_LIBS_BUILD_CMND)
	$(BELOS_ALL_LIBS_BUILD_CMND)
	$(ANASAZI_ALL_LIBS_BUILD_CMND)
	$(RBGEN_ALL_LIBS_BUILD_CMND)
	$(PLIRIS_ALL_LIBS_BUILD_CMND)
	$(ML_ALL_LIBS_BUILD_CMND)
	$(STRATIMIKOS_ALL_LIBS_BUILD_CMND)
	$(MOERTEL_ALL_LIBS_BUILD_CMND)
	$(MEROS_ALL_LIBS_BUILD_CMND)
	$(NOX_ALL_LIBS_BUILD_CMND)
	$(CAPO_ALL_LIBS_BUILD_CMND)
	$(RYTHMOS_ALL_LIBS_BUILD_CMND)
	$(NEW_PACKAGE_ALL_LIBS_BUILD_CMND)
	$(MOOCHO_ALL_LIBS_BUILD_CMND)
	$(SACADO_ALL_LIBS_BUILD_CMND)
	$(PYTRILINOS_ALL_LIBS_BUILD_CMND)
	$(DIDASKO_ALL_LIBS_BUILD_CMND)
	$(WEBTRILINOS_ALL_LIBS_BUILD_CMND)
	$(EXTERNAL_ALL_LIBS_BUILD_CMND)

