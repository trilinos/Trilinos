# ************************************************************************
# 
#                Pliris: Parallel Dense Solver Package
#                 Copyright (2004) Sandia Corporation
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

include Make.inc

force_build:

#
#      Libraries of different precisions for the solver
#       s-single precision real
#       d-double precision real
#       c-complex
#       z- double complex


slib: force_build
	$(MAKE) -f Make.typed lib pre=s typ=SREAL

dlib: force_build
	$(MAKE) -f Make.typed lib pre=d typ=DREAL

clib: force_build
	$(MAKE) -f Make.typed lib pre=c typ=SCPLX

zlib: force_build
	$(MAKE) -f Make.typed lib pre=z typ=ZCPLX

#
#      Test routines to exercise the solvers
#       dftest - fortran double precision
#       stest - c-code single precision real
#       dtest - c-code double precision real
#       ctest - c-code complex
#       ztest - c-code double complex
#

dftest:force_build
	$(MAKE) -f Make.typed dftest pre=d typ=DREAL 

stest: force_build
	$(MAKE) -f Make.typed test  pre=s typ=SREAL

dtest: force_build
	$(MAKE) -f Make.typed test  pre=d typ=DREAL

ctest: force_build
	$(MAKE) -f Make.typed test  pre=c typ=SCPLX

ztest:force_build
	$(MAKE) -f Make.typed test  pre=z typ=ZCPLX

dtest_hard: force_build
	$(MAKE) -f Make.typed hard_test  pre=d typ=DREAL

steve_test: force_build
	$(MAKE) -f Make.typed steve_test  pre=d typ=DREAL
