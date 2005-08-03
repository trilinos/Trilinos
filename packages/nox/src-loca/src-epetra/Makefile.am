## $Id$ 
## $Source$ 
# @HEADER
# ************************************************************************
# 
#                  LOCA Continuation Algorithm Package
#                 Copyright (2005) Sandia Corporation
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
# Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
# (rppawlo@sandia.gov), Sandia National Laboratories.
# 
# ************************************************************************
# @HEADER


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside NOX. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

LOCA_EPETRA_H = \
	$(srcdir)/LOCA_Epetra.H \
	$(srcdir)/LOCA_Epetra_Group.H \
	$(srcdir)/LOCA_Epetra_Interface.H \
	$(srcdir)/LOCA_Epetra_BorderedOp.H \
	$(srcdir)/LOCA_Epetra_HouseholderJacOp.H \
	$(srcdir)/LOCA_Epetra_ShiftInvertOperator.H

LOCA_EPETRA = \
	$(srcdir)/LOCA_Epetra_Group.C \
	$(srcdir)/LOCA_Epetra_BorderedOp.C \
	$(srcdir)/LOCA_Epetra_HouseholderJacOp.C \
	$(srcdir)/LOCA_Epetra_ShiftInvertOperator.C

LOCA_EPETRANEW_H = \
	$(srcdir)/LOCA_EpetraNew_Group.H \
	$(srcdir)/LOCA_EpetraNew_Interface_Required.H \
	$(srcdir)/LOCA_EpetraNew_Interface_xyzt.H \
	$(srcdir)/LOCA_EpetraNew_Interface_MassMatrix.H \
        $(srcdir)/LOCA_EpetraNew_Interface_TimeDependent.H \
	$(srcdir)/LOCA_Epetra_CompactWYOp.H \
	$(srcdir)/LOCA_BorderedSystem_EpetraHouseholder.H \
	$(srcdir)/LOCA_Epetra_AugmentedOp.H \
	$(srcdir)/LOCA_BorderedSystem_EpetraAugmented.H \
	$(srcdir)/LOCA_Epetra_Factory.H
 
LOCA_EPETRANEW = \
	$(srcdir)/LOCA_EpetraNew_Group.C \
	$(srcdir)/LOCA_EpetraNew_Interface_xyzt.C \
	$(srcdir)/LOCA_Epetra_CompactWYOp.C \
	$(srcdir)/LOCA_BorderedSystem_EpetraHouseholder.C \
	$(srcdir)/LOCA_Epetra_AugmentedOp.C \
	$(srcdir)/LOCA_BorderedSystem_EpetraAugmented.C \
	$(srcdir)/LOCA_Epetra_Factory.C

# ------------------------------------------------------------------------
# Includes for required and optional library dependencies
# ------------------------------------------------------------------------
include $(top_builddir)/Makefile.export.nox

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(LOCA_INCLUDES))
else
EXPORT_INCLUDES = $(LOCA_INCLUDES)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

# ------------------------------------------------------------------------
# Noxepetra library specifications
# ------------------------------------------------------------------------

liblocaepetra_la_SOURCES = \
	$(LOCA_EPETRA) \
	$(LOCA_EPETRA_H) \
	$(LOCA_EPETRANEW) \
	$(LOCA_EPETRANEW_H)

#liblocaepetra_la_CXXFLAGS = \
#	-I$(top_srcdir)/src-loca/src \
#	-I$(top_srcdir)/src \
#	-I$(top_srcdir)/src-epetra \
#	$(TRILINOS_CXXFLAGS) $(LIBTOOL_COMPILER_FLAGS)

#liblocaepetra_la_LIBADD = \
#	$(top_builddir)/src-loca/src/libloca.la \
#	$(CXXREPOS) $(XTRALDADD) $(LIBTOOL_CCLINKER_FLAGS)

# ------------------------------------------------------------------------
# When to turn it on
# ------------------------------------------------------------------------

if BUILD_LIBLOCAEPETRA

# List of all libraries to install in $(libexecdir)

lib_LTLIBRARIES = liblocaepetra.la

# List of all headers to install in $(includedir)

include_HEADERS = $(LOCA_EPETRA_H) $(LOCA_EPETRANEW_H)

endif

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

# ------------------------------------------------------------------------
# Files to be deleted by 'make clean'
# ------------------------------------------------------------------------

MOSTLYCLEANFILES = $(CXXREPOS)

# ------------------------------------------------------------------------
# Include global libtool definitions
# ------------------------------------------------------------------------
include $(top_srcdir)/config/libtool_macros.am




