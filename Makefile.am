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
EXTRA_DIST = \
	testharness/subdir-list testharness/checkin-test-harness \
	sampleScripts/aix_mpi_xlC sampleScripts/aix_mpi_xlC2 \
	sampleScripts/compaqCluster_mpi_stratus \
	sampleScripts/compaqCluster_serial_stratus sampleScripts/cplant_mpi_ross \
	sampleScripts/cplant_serial_ross \
	sampleScripts/cygwin_mpi sampleScripts/cygwin_serial_gcc \
	sampleScripts/cygwin_serial_icl sampleScripts/cygwingcc_amesos \
	sampleScripts/irix_mpi_sierra sampleScripts/janus_mpi \
	sampleScripts/janus_mpi_stlport \
	sampleScripts/linux_mpi_codecov \
	sampleScripts/linux_mpi_herouxsmp sampleScripts/linux_mpi_instCluster \
	sampleScripts/linux_serial_codecov \
	sampleScripts/macosx_mpi_paul sampleScripts/macosx_serial_paul \
	sampleScripts/macosx_serial_samt5980 \
	sampleScripts/scyld_mpi_beowulf sampleScripts/sgi64_mpi_atlantis \
	sampleScripts/sgi64_serial_atlantis \
	sampleScripts/solaris_mpi_conejo sampleScripts/solaris_serial_paul \
	sampleScripts/README config/AutotoolsVersionInfo

## #######################################################################
## Tools in the auxillary directory 
## #######################################################################
AUX_DIST = config/install-sh config/missing config/mkinstalldirs

## #######################################################################
## Files to be deleted by 'make maintainer-clean'
## #######################################################################
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 autom4te.cache/* \
	configure config.status config.log \
	src/common/config-h.in src/common/stamp-h.in \
	$(AUX_DIST) 

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################
SUBDIRS = packages
