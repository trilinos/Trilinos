## $Id$ 
## $Source$ 
# @HEADER
# ************************************************************************
# 
#            NOX: An Object-Oriented Nonlinear Solver Package
#                 Copyright (2002) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#   
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#   
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
# 
# Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
# (rppawlo@sandia.gov).
# 
# ************************************************************************
# @HEADER


# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside NOX. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

LOCA_ABSTRACT_H = \
        $(srcdir)/LOCA_Abstract_Group.H \
	$(srcdir)/LOCA_Abstract_DataOutput.H

LOCA_ABSTRACT = \
        $(srcdir)/LOCA_Abstract_Group.C

LOCA_UTILS_H = \
        $(srcdir)/LOCA_Utils.H

LOCA_UTILS = \
        $(srcdir)/LOCA_Utils.C

LOCA_PARAMETER_VECTOR_H = \
        $(srcdir)/LOCA_Parameter_Vector.H

LOCA_PARAMETER_VECTOR = \
        $(srcdir)/LOCA_Parameter_Vector.C

LOCA_STEPPER_H = \
	$(srcdir)/LOCA_Stepper.H \
	$(srcdir)/LOCA_StepperArclength.H

LOCA_STEPPER = \
	$(srcdir)/LOCA_Stepper.C \
	$(srcdir)/LOCA_StepperArclength.C

LOCA_BIFURCATION_H = \
	$(srcdir)/LOCA_MultiVector.H \
        $(srcdir)/LOCA_Bifurcation_TPBordVector.H \
        $(srcdir)/LOCA_Bifurcation_TPBordGroup.H \
        $(srcdir)/LOCA_Bifurcation_ArcLengthVector.H \
        $(srcdir)/LOCA_Bifurcation_ArcLengthGroup.H

LOCA_BIFURCATION = \
	$(srcdir)/LOCA_MultiVector.C \
        $(srcdir)/LOCA_Bifurcation_TPBordVector.C \
        $(srcdir)/LOCA_Bifurcation_TPBordGroup.C \
        $(srcdir)/LOCA_Bifurcation_ArcLengthVector.C \
        $(srcdir)/LOCA_Bifurcation_ArcLengthGroup.C

LOCA_DERIVUTILS_H = \
        $(srcdir)/LOCA_DerivUtils.H

LOCA_DERIVUTILS = \
        $(srcdir)/LOCA_DerivUtils.C

LOCA_H = \
        $(LOCA_ABSTRACT_H) \
	$(LOCA_UTILS_H) \
	$(LOCA_PARAMETER_VECTOR_H) \
	$(LOCA_STEPPER_H) \
        $(LOCA_SOLVER_H) \
	$(LOCA_BIFURCATION_H) \
	$(LOCA_DERIVUTILS_H)

LOCA = \
	$(LOCA_ABSTRACT) \
	$(LOCA_UTILS) \
	$(LOCA_PARAMETER_VECTOR) \
	$(LOCA_STEPPER) \
        $(LOCA_SOLVER) \
	$(LOCA_BIFURCATION) \
	$(LOCA_DERIVUTILS)

# ------------------------------------------------------------------------
# Loca library specifications
# ------------------------------------------------------------------------

libloca_a_SOURCES = $(LOCA_H) $(LOCA)

libloca_a_CXXFLAGS = -I$(top_srcdir)/src

# ------------------------------------------------------------------------
# When to turn it on
# ------------------------------------------------------------------------

if BUILD_LIBLOCA

# List of all libraries to install in $(libexecdir)

lib_LIBRARIES = libloca.a

# List of all headers to install in $(pkgincludedir)

pkginclude_HEADERS = $(LOCA_H)

endif

# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libloca_a_AR = $(ALTERNATE_AR)

else

libloca_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libloca_a_LIBADD = $(CXXREPOS)

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

