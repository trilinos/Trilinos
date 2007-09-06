#@HEADER
# ************************************************************************
# 
#        AztecOO: An Object-Oriented Aztec Linear Solver Package 
#                 Copyright (2002) Sandia Corporation
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
#@HEADER

# The following line helps the test harness recover from build errors.

all-local:

#SUBDIRS = azoo_iterate_hb aztec_app aztec_hb AztecOO AztecOO_UserOpUserMat \
#	IfpackAztecOO simple_epetra_aztecoo
# note that simple_epetra_aztecoo, any AztecOO*, and azoo should not be built when "core Aztec" is 
# being built

#if BUILD_TRIUTILS
USER_OP_USER_MAT_SUBDIR=AztecOO_UserOpUserMat
AZTECOO_SUBDIR=AztecOO
#endif

if HAVE_AZTECOO_EPETRAEXT
AZTECOO_EPETRAEXT_SUBDIR=AztecOO_MatlabInput
endif

SUBDIRS = $(AZTECOO_SUBDIR) $(USER_OP_USER_MAT_SUBDIR) $(AZTECOO_EPETRAEXT_SUBDIR) simple_epetra_aztecoo UG_Ex1 ClassicAztec
