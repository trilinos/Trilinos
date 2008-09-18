## $Id$ 
## $Source$ 
#@HEADER
# ************************************************************************
# 
#            NOX: An Object-Oriented Nonlinear Solver Package
#                 Copyright (2002) Sandia Corporation
# 
#            LOCA: Library of Continuation Algorithms Package
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
# 
# Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
# Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
# ************************************************************************
#@HEADER

# ------------------------------------------------------------------------
# FEApp Sources
# ------------------------------------------------------------------------
include $(top_srcdir)/../sacado/example/FEApp/MakefileFEAppStub.am

# Include export Makefiles for correct paths
include $(top_builddir)/Makefile.export.nox

EXEEXT = .exe

if BUILD_EPETRA_EXAMPLES
if ENABLE_LOCA_SACADO
if BUILD_NOX_EPETRAEXT_SUPPORT
include $(top_builddir)/../sacado/Makefile.export.sacado
noinst_LIBRARIES = libfeapp.a
noinst_PROGRAMS = brusselator
endif
endif
endif

AM_CPPFLAGS = -I$(FEAPP_PATH) $(NOX_INCLUDES) $(SACADO_INCLUDES)

libfeapp_a_SOURCES = $(FEAPP_H) $(FEAPP)

if USE_ALTERNATE_AR
libfeapp_a_AR = $(ALTERNATE_AR)
else
libfeapp_a_AR = $(AR) cru
endif

brusselator_SOURCES = $(srcdir)/brusselator.C 

brusselator_DEPENDENCIES = $(builddir)/libfeapp.a

brusselator_LDADD = -L$(builddir) -lfeapp $(NOX_LIBS) $(SACADO_LIBS)