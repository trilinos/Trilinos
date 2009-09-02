#@HEADER
# ***********************************************************************
#
#                           Rythmos Package
#                 Copyright (2006) Sandia Corporation
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
# Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
#
# ***********************************************************************
#@HEADER

# The following line helps the test harness recover from build errors.

all-local:

if ENABLE_RYTHMOS_STRATIMIKOS
if ENABLE_RYTHMOS_THYRA_EPETRAEXT
RYTHMOS_STRATIMIKOS_EXAMPLES = diagonalTransient timeDiscretizedBackwardEuler
if ENABLE_RYTHMOS_NOX
if ENABLE_RYTHMOS_NOX_EPETRA
RYTHMOS_NOX_EXAMPLES = epetra
endif # ENABLE_RYTHMOS_NOX_EPETRA
endif # ENABLE_RYTHMOS_NOX
endif # ENABLE_RYTHMOS_THYRA_EPETRAEXT
endif # ENABLE_RYTHMOS_STRATIMIKOS

SUBDIRS = $(RYTHMOS_STRATIMIKOS_EXAMPLES) $(RYTHMOS_NOX_EXAMPLES) 

