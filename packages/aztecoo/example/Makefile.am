#@HEADER
# ************************************************************************
# 
#            Trilinos: An Object-Oriented Solver Framework
#                Copyright (2002) Sandia Corporation
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
# 
# ************************************************************************
#@HEADER

#SUBDIRS = azoo_iterate_hb aztec_app aztec_hb AztecOO AztecOO_UserOpUserMat \
#	IfpackAztecOO simple_epetra_aztecoo
# note that simple_epetra_aztecoo, any AztecOO*, and azoo should not be built when "core Aztec" is 
# being built

#if BUILD_TRIUTILS
USER_OP_USER_MAT_SUBDIR=AztecOO_UserOpUserMat
AZTECOO_SUBDIR=AztecOO
#endif

SUBDIRS = $(AZTECOO_SUBDIR) $(USER_OP_USER_MAT_SUBDIR) simple_epetra_aztecoo
