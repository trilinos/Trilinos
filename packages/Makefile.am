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
ACLOCAL_AMFLAGS = -I config

## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
EXTRA_DIST = common/Doxyfile common/header.html common/footer.html \
common/styles.css common/tag_files/keepme.txt CMakeLists.txt

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
## Names of package subdirectories to be make'd recursively
## #######################################################################

# 2007/07/18: rabartl: Note: Adding or removing a package now only requires
# changing three lines!

if BUILD_AMESOS
AMESOS_SUBDIR = amesos
endif

if BUILD_ANASAZI
ANASAZI_SUBDIR = anasazi
endif

if BUILD_AZTECOO
AZTECOO_SUBDIR = aztecoo
endif

if BUILD_BELOS
BELOS_SUBDIR = belos
endif

if BUILD_CLAPS
CLAPS_SUBDIR = claps
endif

if BUILD_DIDASKO
DIDASKO_SUBDIR = didasko
endif

if BUILD_EPETRA
EPETRA_SUBDIR = epetra
endif

if BUILD_EPETRAEXT
EPETRAEXT_SUBDIR = epetraext
endif

if BUILD_EXTERNAL
EXTERNAL_SUBDIR = external
endif

if BUILD_FEI
FEI_SUBDIR = fei
endif

if BUILD_IFPACK
IFPACK_SUBDIR = ifpack
endif

if BUILD_INTREPID
INTREPID_SUBDIR = intrepid
endif

if BUILD_GALERI
GALERI_SUBDIR = galeri 
endif

if BUILD_KOKKOS
KOKKOS_SUBDIR = kokkos
endif

if BUILD_KOMPLEX
KOMPLEX_SUBDIR = komplex
endif

if BUILD_MEROS
MEROS_SUBDIR = meros
endif

if BUILD_ML
ML_SUBDIR = ml
endif

if BUILD_MOERTEL
MOERTEL_SUBDIR = moertel 
endif

if BUILD_MOOCHO
MOOCHO_SUBDIR = moocho
endif

if BUILD_PAMGEN
PAMGEN_SUBDIR = pamgen
endif

if BUILD_NEW_PACKAGE
NEW_PACKAGE_SUBDIR = new_package
endif

if BUILD_ISORROPIA
ISORROPIA_SUBDIR = isorropia
endif

if BUILD_NOX
NOX_SUBDIR = nox
endif

if BUILD_PHALANX
PHALANX_SUBDIR = phalanx
endif

if BUILD_PHDMESH
PHDMESH_SUBDIR = phdmesh
endif

if BUILD_PLIRIS
PLIRIS_SUBDIR = pliris
endif

if BUILD_PYTRILINOS
PYTRILINOS_SUBDIR = PyTrilinos 
endif

if BUILD_RBGEN
RBGEN_SUBDIR = rbgen
endif

if BUILD_RYTHMOS
RYTHMOS_SUBDIR = rythmos
endif

if BUILD_SACADO
SACADO_SUBDIR = sacado
endif

if BUILD_TEUCHOS
TEUCHOS_SUBDIR = teuchos
endif

if BUILD_RTOP
RTOP_SUBDIR = rtop
endif

if BUILD_THREADPOOL
THREADPOOL_SUBDIR = ThreadPool
endif

if BUILD_THYRA
THYRA_SUBDIR = thyra
endif

if BUILD_STRATIMIKOS
STRATIMIKOS_SUBDIR = stratimikos
endif

if BUILD_SUNDANCE
SUNDANCE_SUBDIR = Sundance
endif

if BUILD_TPETRA
TPETRA_SUBDIR = tpetra
endif

if BUILD_TRILINOSCOUPLINGS
TRILINOSCOUPLINGS_SUBDIR = trilinoscouplings
endif

if BUILD_TRIUTILS
TRIUTILS_SUBDIR = triutils
endif

if BUILD_WEBTRILINOS
WEBTRILINOS_SUBDIR = WebTrilinos
endif

if BUILD_ZOLTAN
ZOLTAN_SUBDIR = zoltan
endif

## #######################################################################
## Single list of configured pacakges
## #######################################################################

#
# This list of pacakges to be built are roughly in the order of lower-level to
# higher-level pacakges.
#
# Listing the packages from lower level to higher level allows the test
# harness to perform more efficiently as it disables broken packages and
# therefore disables higher-level pacakges.  In general, these packages should
# be listed in the reverse order that the libraries would be listed on the
# link line.
#

PACKAGE_NAMES = \
	$(TEUCHOS_SUBDIR)      \
	$(RTOP_SUBDIR)         \
	$(KOKKOS_SUBDIR)       \
	$(ZOLTAN_SUBDIR)       \
	$(EPETRA_SUBDIR)       \
	$(TRIUTILS_SUBDIR)     \
	$(TPETRA_SUBDIR)       \
	$(EPETRAEXT_SUBDIR)    \
	$(THYRA_SUBDIR)        \
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
	$(INTREPID_SUBDIR)     \
	$(MOERTEL_SUBDIR)      \
	$(STRATIMIKOS_SUBDIR)  \
	$(MEROS_SUBDIR)        \
	$(FEI_SUBDIR)        \
	$(SACADO_SUBDIR)       \
	$(NOX_SUBDIR)          \
	$(PHALANX_SUBDIR)      \
	$(RYTHMOS_SUBDIR)      \
	$(TRILINOSCOUPLINGS_SUBDIR) \
	$(NEW_PACKAGE_SUBDIR)  \
	$(MOOCHO_SUBDIR)       \
	$(PAMGEN_SUBDIR)      \
	$(THREADPOOL_SUBDIR) \
	$(PHDMESH_SUBDIR) \
	$(SUNDANCE_SUBDIR)     \
	$(PYTRILINOS_SUBDIR)   \
	$(DIDASKO_SUBDIR) \
	$(WEBTRILINOS_SUBDIR) \
	$(EXTERNAL_SUBDIR)

SUBDIRS = $(PACKAGE_NAMES)


#
# Standard build targets which build tests and examples
#

examples:
	$(top_srcdir)/config/makeall $(MAKE) $@ $(PACKAGE_NAMES)

install-examples:
	$(top_srcdir)/config/makeall $(MAKE) $@ $(PACKAGE_NAMES)

clean-examples:
	$(top_srcdir)/config/makeall $(MAKE) $@ $(PACKAGE_NAMES)

tests:
	$(top_srcdir)/config/makeall $(MAKE) $@ $(PACKAGE_NAMES)

clean-tests:
	$(top_srcdir)/config/makeall $(MAKE) $@ $(PACKAGE_NAMES)

everything:
	 $(MAKE) && $(MAKE) examples && $(MAKE) tests

install-everything:
	$(MAKE) install && $(MAKE) install-examples && $(MAKE) install-finish

clean-everything:
	$(MAKE) clean-examples && $(MAKE) clean-tests && $(MAKE) clean

#The below target was added because make install was previously not called in
#the top level directory when using make install-everything, which caused 
#these files to be missing from the install directory.
install-finish:
	cp $(top_srcdir)/Trilinos_version.h $(DESTDIR)$(includedir)/. ;
	cp $(top_srcdir)/commonTools/buildTools/strip_dup_incl_paths.pl $(DESTDIR)$(includedir)/. ;
	cp $(top_srcdir)/commonTools/buildTools/strip_dup_libs.pl $(DESTDIR)$(includedir)/.

partial-clean:
	find . -name "*.a" -exec rm {} \;
	find . -name "*.exe" -exec rm {} \; 
	find . -name "*.xml" -exec rm {} \; 
	find . -name "*.pl" -exec rm {} \; 
	find . -name "*.csh" -exec rm {} \; 
	find . -name "*.sh" -exec rm {} \; 
	find . -name "*.cua" -exec rm {} \; 
	find . -name "*.rua" -exec rm {} \; 
	find . -name "*.hb" -exec rm {} \; 
	find . -name "*.mtx" -exec rm {} \; 
	find . -name "*.mm" -exec rm {} \; 
	find . -name "*.triU" -exec rm {} \; 
	find . -name "*.triS" -exec rm {} \; 
	find . -name "*.dat" -exec rm {} \; 
	find . -name "*.m" -exec rm {} \; 
	find . -name "*.grid" -exec rm {} \; 
	find . -name "*.txt" -exec rm {} \; 
	find . -name "*.evaluateCriteria" -exec rm {} \; 

# 2007/07/18: rabartl: Above, a simple script makeall has to be used since a
# for loop with shell variables is employed.  I tried to put the for look
# directly in the make target but 'make' simply replaced the shell varibles
# names (with empty values) before the shell command was even executed.  This
# made it necessary to call a shell script so that the shell variables would
# not be replaced.  The use of the script 'makeall' also results in very clean
# target rules as shown above and improves maintenance.
